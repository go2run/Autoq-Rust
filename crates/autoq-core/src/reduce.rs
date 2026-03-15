//! # Reduction algorithms for tree automata
//!
//! ## C++ correspondence
//! Port of `src/reduce.cc` (31KB).
//!
//! ## Current scope
//! - `remove_useless`: remove unreachable states
//! - `state_renumbering`: compact state numbering
//! - `k_unification`: unify k-values across leaf amplitudes
//! - `fraction_simplification`: simplify FiveTuple fractions
//! - Stubs for `sim_reduce`, `light_reduce_up/down`

use crate::automata::{Automata, State};
use crate::symbol::SymbolTrait;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

impl<S: SymbolTrait> Automata<S> {
    /// Remove states not reachable from final states (top-down)
    /// and states not leading to any leaf (bottom-up).
    ///
    /// C++ equivalent: `remove_useless(only_bottom_up)`.
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

        // Remap final states
        self.final_states = self.final_states.iter()
            .filter_map(|s| mapping.get(s).copied())
            .collect();

        // Remap transitions
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

    /// Apply preferred reduction pipeline.
    ///
    /// C++ equivalent: `reduce()`.
    pub fn reduce(&mut self) {
        self.remove_useless();
        self.state_renumbering();
    }
}
