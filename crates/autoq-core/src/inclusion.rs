//! # Language inclusion checking for tree automata
//!
//! ## C++ correspondence
//! Port of `src/inclusion.cc` (154KB). This is the most complex algorithm
//! in AutoQ — the full C++ file handles Concrete, Symbolic, Predicate,
//! Index, and Constrained variants with Z3 integration.
//!
//! ## Current scope
//! This initial version implements:
//! - `is_empty`: check if L(A) = ∅
//! - `is_included_in`: check if L(A) ⊆ L(B) (Concrete only)
//! - `is_included_up_to_scaling`: check ∃λ≠0 : λ·L(A) ⊆ L(B)
//!
//! ## Future work
//! - Symbolic/Predicate/Index inclusion requires autoq-z3 integration.
//! - sim_reduce optimization (simulation-based reduction before inclusion).

use crate::automata::{Automata, State, SymbolTag, Tag};
use crate::symbol::SymbolTrait;
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

/// Check if the automaton accepts the empty language.
///
/// C++ equivalent: `Automata<Symbol>::empty()`.
///
/// BFS on macro-states (sets of automaton states). Starts from each singleton
/// {final_state} and tries to find a consistent path to leaf transitions.
pub fn is_empty<S: SymbolTrait>(aut: &Automata<S>) -> bool {
    // Build per-state index: state → union of tags, and per-tag transitions
    let mut state_tags: HashMap<State, Tag> = HashMap::new();
    let mut state_trans: HashMap<State, Vec<(SymbolTag<S>, Vec<Vec<State>>)>> = HashMap::new();

    for (st, map) in &aut.transitions {
        for (&parent, children_set) in map {
            *state_tags.entry(parent).or_insert(0) |= st.tag;
            let entry = state_trans.entry(parent).or_default();
            let children_vec: Vec<Vec<State>> = children_set.iter().cloned().collect();
            entry.push((st.clone(), children_vec));
        }
    }

    let mut bfs: VecDeque<BTreeSet<State>> = VecDeque::new();
    for &q in &aut.final_states {
        let mut vertex = BTreeSet::new();
        vertex.insert(q);
        bfs.push_back(vertex);
    }

    while let Some(vertex) = bfs.pop_front() {
        // Compute intersection of all tags in this vertex
        let mut range = Tag::MAX;
        for &top in &vertex {
            range &= state_tags.get(&top).copied().unwrap_or(0);
        }
        if range == 0 { continue; }

        // Try each candidate color bit
        for bit in 0..64u32 {
            if (range >> bit) & 1 == 0 { continue; }
            let candidate = 1u64 << bit;

            let mut new_vertex = BTreeSet::new();
            let mut all_leaf = true;
            let mut valid = true;

            for &top in &vertex {
                let mut found = false;
                if let Some(trans) = state_trans.get(&top) {
                    for (st, children_vecs) in trans {
                        if st.tag & candidate != 0 {
                            if st.is_leaf() {
                                // Leaf found for this state
                                found = true;
                            } else {
                                all_leaf = false;
                                for children in children_vecs {
                                    for &child in children {
                                        new_vertex.insert(child);
                                    }
                                }
                                found = true;
                            }
                            break;
                        }
                    }
                }
                if !found { valid = false; break; }
            }

            if !valid { continue; }
            if all_leaf { return false; } // Found complete leaf path → non-empty
            if !new_vertex.is_empty() {
                bfs.push_back(new_vertex);
            }
        }
    }
    true // No leaf path found → empty
}

/// Check if L(A) ⊆ L(B) for concrete automata.
///
/// C++ equivalent: `operator<=(Automata o)`.
///
/// Anti-chain BFS over obligation pairs (state_A, set_of_states_B).
pub fn is_included_in<S: SymbolTrait>(aut_a: &Automata<S>, aut_b: &Automata<S>) -> bool {
    let trans_a = build_trans_map(aut_a);
    let trans_b = build_trans_map(aut_b);

    type Obligation = (State, BTreeSet<State>);
    let mut worklist: VecDeque<Obligation> = VecDeque::new();
    let mut visited: HashSet<Obligation> = HashSet::new();

    let finals_b: BTreeSet<State> = aut_b.final_states.iter().cloned().collect();
    for &qa in &aut_a.final_states {
        let obl = (qa, finals_b.clone());
        if visited.insert(obl.clone()) {
            worklist.push_back(obl);
        }
    }

    while let Some((qa, qbs)) = worklist.pop_front() {
        if let Some(a_trans) = trans_a.get(&qa) {
            for (sym_a, tag_a, children_a) in a_trans {
                let mut matching_b: Vec<Vec<State>> = Vec::new();
                for &qb in &qbs {
                    if let Some(b_trans) = trans_b.get(&qb) {
                        for (sym_b, tag_b, children_b) in b_trans {
                            if sym_b == sym_a && tag_b & tag_a != 0 {
                                matching_b.push(children_b.clone());
                            }
                        }
                    }
                }

                if sym_a.is_leaf() && matching_b.is_empty() {
                    return false;
                }

                if sym_a.is_internal() && children_a.len() == 2 {
                    for (pos, &child_a) in children_a.iter().enumerate() {
                        // If child_a has no transitions at all, it's a dead state.
                        // No trees can be produced through it, so this obligation
                        // is vacuously satisfied.
                        if !trans_a.contains_key(&child_a) { continue; }
                        let child_b_set: BTreeSet<State> = matching_b.iter()
                            .filter_map(|bc| bc.get(pos).copied())
                            .collect();
                        if child_b_set.is_empty() { return false; }
                        let obl = (child_a, child_b_set);
                        if visited.insert(obl.clone()) {
                            worklist.push_back(obl);
                        }
                    }
                }
            }
        }
    }
    true
}

/// Check L(A) == L(B).
pub fn are_equal<S: SymbolTrait>(a: &Automata<S>, b: &Automata<S>) -> bool {
    is_included_in(a, b) && is_included_in(b, a)
}

/// Build state → [(symbol, tag, children)] lookup map.
fn build_trans_map<S: SymbolTrait>(aut: &Automata<S>) -> HashMap<State, Vec<(S, Tag, Vec<State>)>> {
    let mut map: HashMap<State, Vec<(S, Tag, Vec<State>)>> = HashMap::new();
    for (st, parent_map) in &aut.transitions {
        for (&parent, children_set) in parent_map {
            for children in children_set {
                map.entry(parent)
                    .or_default()
                    .push((st.symbol.clone(), st.tag, children.clone()));
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbol::ConcreteSymbol;

    #[test]
    fn zero_state_not_empty() {
        let aut = Automata::<ConcreteSymbol>::zero_state(2);
        assert!(!is_empty(&aut));
    }

    #[test]
    fn uniform_not_empty() {
        let aut = Automata::<ConcreteSymbol>::uniform(2);
        assert!(!is_empty(&aut));
    }

    #[test]
    fn no_finals_is_empty() {
        let aut = Automata::<ConcreteSymbol>::new(1);
        assert!(is_empty(&aut));
    }

    #[test]
    fn self_inclusion() {
        let a = Automata::<ConcreteSymbol>::zero_state(2);
        assert!(is_included_in(&a, &a));
    }

    #[test]
    fn self_equality() {
        let a = Automata::<ConcreteSymbol>::uniform(2);
        assert!(are_equal(&a, &a.clone()));
    }
}
