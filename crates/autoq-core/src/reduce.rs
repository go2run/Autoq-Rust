//! # Reduction algorithms for tree automata
//!
//! ## C++ correspondence
//! Port of `src/reduce.cc` (31KB).
//!
//! ## Algorithms
//! - `remove_useless`: remove states unreachable from final states (top-down)
//! - `state_renumbering`: compact state numbering
//! - `bottom_up_reduce`: bottom-up equivalence class merging (the main reducer)

use crate::automata::{Automata, State, SymbolTag, Tag, TopDownTransitions};
use crate::symbol::SymbolTrait;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

impl<S: SymbolTrait> Automata<S> {
    /// Remove states not reachable from final states (top-down).
    ///
    /// C++ equivalent: `remove_useless()`.
    pub fn remove_useless(&mut self) {
        // Top-down: find states reachable from final states
        let mut reachable = BTreeSet::new();
        let mut queue: VecDeque<State> = self.final_states.iter().cloned().collect();
        while let Some(s) = queue.pop_front() {
            if reachable.insert(s) {
                for (_, map) in &self.transitions {
                    if let Some(children_set) = map.get(&s) {
                        for children in children_set {
                            for &child in children {
                                queue.push_back(child);
                            }
                        }
                    }
                }
            }
        }

        // Remove transitions with unreachable parent states
        self.transitions.retain(|_, map| {
            map.retain(|parent, _| reachable.contains(parent));
            !map.is_empty()
        });
    }

    /// Compact state numbering: reassign state IDs to 0..n-1.
    ///
    /// C++ equivalent: `state_renumbering()`.
    pub fn state_renumbering(&mut self) {
        let mut all_states = BTreeSet::new();
        for &s in &self.final_states { all_states.insert(s); }
        for (_, map) in &self.transitions {
            for (&parent, children_set) in map {
                all_states.insert(parent);
                for children in children_set {
                    for &child in children { all_states.insert(child); }
                }
            }
        }

        let mapping: BTreeMap<State, State> = all_states.iter()
            .enumerate()
            .map(|(new, &old)| (old, new as State))
            .collect();

        self.final_states = self.final_states.iter()
            .filter_map(|s| mapping.get(s).copied())
            .collect();

        let old_trans = std::mem::take(&mut self.transitions);
        for (st, map) in old_trans {
            for (parent, children_set) in map {
                let new_parent = mapping[&parent];
                for children in children_set {
                    let new_children: Vec<State> = children.iter()
                        .map(|c| mapping[c])
                        .collect();
                    self.transitions.entry(st.clone()).or_default()
                        .entry(new_parent).or_default()
                        .insert(new_children);
                }
            }
        }

        self.state_num = all_states.len() as State;
    }

    /// Bottom-up equivalence class reduction.
    ///
    /// C++ equivalent: `bottom_up_reduce()`.
    ///
    /// Iterates from leaves upward, layer by layer. At each layer,
    /// states with identical "downward signatures" (same transitions below them)
    /// are merged into a single representative state.
    pub fn bottom_up_reduce(&mut self) {
        // Signature: for each state, the set of (Symbol → (children → tag))
        type StateSig<S> = BTreeMap<S, BTreeMap<Vec<State>, Tag>>;

        let mut transitions2: TopDownTransitions<S> = BTreeMap::new();
        let mut rewrite: BTreeMap<State, State> = BTreeMap::new();

        let resolve = |s: State, rw: &BTreeMap<State, State>| -> State {
            *rw.get(&s).unwrap_or(&s)
        };

        // Separate leaf and internal transitions
        let mut leaf_trans: Vec<(SymbolTag<S>, BTreeMap<State, BTreeSet<Vec<State>>>)> = Vec::new();
        let mut internal_by_qubit: BTreeMap<i64, Vec<(SymbolTag<S>, BTreeMap<State, BTreeSet<Vec<State>>>)>> = BTreeMap::new();

        for (st, map) in &self.transitions {
            if st.is_leaf() {
                leaf_trans.push((st.clone(), map.clone()));
            } else {
                internal_by_qubit.entry(st.symbol.qubit())
                    .or_default()
                    .push((st.clone(), map.clone()));
            }
        }

        // Process leaves: build signatures
        {
            let mut qfic: BTreeMap<State, StateSig<S>> = BTreeMap::new();
            for (st, map) in &leaf_trans {
                for (&parent, children_set) in map {
                    for children in children_set {
                        *qfic.entry(parent).or_default()
                            .entry(st.symbol.clone()).or_default()
                            .entry(children.clone()).or_default() |= st.tag;
                    }
                }
            }

            let mut repr: BTreeMap<StateSig<S>, State> = BTreeMap::new();
            for (top, sig) in &qfic {
                if let Some(&canonical) = repr.get(sig) {
                    rewrite.insert(*top, canonical);
                } else {
                    rewrite.insert(*top, *top);
                    repr.insert(sig.clone(), *top);
                    for (sym, children_tags) in sig {
                        for (children, &tag) in children_tags {
                            transitions2.entry(SymbolTag::new(sym.clone(), tag)).or_default()
                                .entry(*top).or_default()
                                .insert(children.clone());
                        }
                    }
                }
            }
        }

        // Process internal transitions from highest qubit to lowest
        let qubits: Vec<i64> = internal_by_qubit.keys().rev().cloned().collect();
        for qubit in qubits {
            let layer = internal_by_qubit.get(&qubit).unwrap();

            let mut qfic: BTreeMap<State, StateSig<S>> = BTreeMap::new();
            for (st, map) in layer {
                for (&parent, children_set) in map {
                    for children in children_set {
                        let rewritten: Vec<State> = children.iter()
                            .map(|&s| resolve(s, &rewrite))
                            .collect();
                        *qfic.entry(parent).or_default()
                            .entry(st.symbol.clone()).or_default()
                            .entry(rewritten).or_default() |= st.tag;
                    }
                }
            }

            let mut repr: BTreeMap<StateSig<S>, State> = BTreeMap::new();
            for (top, sig) in &qfic {
                if let Some(&canonical) = repr.get(sig) {
                    rewrite.insert(*top, canonical);
                } else {
                    rewrite.insert(*top, *top);
                    repr.insert(sig.clone(), *top);
                    for (sym, children_tags) in sig {
                        for (children, &tag) in children_tags {
                            transitions2.entry(SymbolTag::new(sym.clone(), tag)).or_default()
                                .entry(*top).or_default()
                                .insert(children.clone());
                        }
                    }
                }
            }
        }

        self.transitions = transitions2;

        // Rewrite final states
        self.final_states = self.final_states.iter()
            .map(|&s| resolve(s, &rewrite))
            .collect();
        self.final_states.sort();
        self.final_states.dedup();

        self.state_renumbering();
    }

    /// Simplify all leaf symbol fractions.
    ///
    /// C++ equivalent: `fraction_simplification()`.
    ///
    /// Ensures equivalent amplitudes with different representations
    /// (e.g., [2,0,0,0,2] and [1,0,0,0,0]) are normalized to the same form.
    pub fn fraction_simplification(&mut self) {
        let old_trans = std::mem::take(&mut self.transitions);
        for (mut st, map) in old_trans {
            if st.is_leaf() {
                S::complex_fraction_simplification(st.symbol.complex_mut());
            }
            self.transitions.entry(st).or_default().extend(map);
        }
    }

    /// Apply preferred reduction pipeline.
    ///
    /// C++ equivalent: `reduce()`.
    pub fn reduce(&mut self) {
        self.fraction_simplification();
        if !self.has_loop {
            self.bottom_up_reduce();
        } else {
            self.remove_useless();
            self.state_renumbering();
        }
    }
}
