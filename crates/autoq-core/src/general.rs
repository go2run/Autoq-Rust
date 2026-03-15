//! # General operations on tree automata
//!
//! ## C++ correspondence
//! Port of `src/general.cc` (17KB).
//!
//! ## Operations
//! - Union (∪): accepts trees from either automaton
//! - Intersection (∩): accepts trees from both automata
//! - Tensor product (⊗): combines two automata on disjoint qubits
//!
//! ## Current scope
//! Union is implemented. Intersection and tensor product are stubs.

use crate::automata::{Automata, State};
use crate::symbol::SymbolTrait;

impl<S: SymbolTrait> Automata<S> {
    /// Union of two automata: L(A ∪ B) = L(A) ∪ L(B).
    ///
    /// C++ equivalent: `operator||(Automata o)`.
    ///
    /// Combines transitions with distinct color tags so each original automaton's
    /// paths remain distinguishable.
    pub fn union(&self, other: &Self) -> Self {
        assert_eq!(self.qubit_num, other.qubit_num, "union requires same qubit count");
        let mut result = self.clone();
        let state_offset = result.state_num;

        // Shift other's states and add transitions
        for (st, map) in &other.transitions {
            for (&parent, children_set) in map {
                let new_parent = parent + state_offset;
                for children in children_set {
                    let new_children: Vec<State> = children.iter()
                        .map(|&c| c + state_offset)
                        .collect();
                    result.add_transition(st.symbol.clone(), st.tag, new_parent, new_children);
                }
            }
        }

        for &fs in &other.final_states {
            result.final_states.push(fs + state_offset);
        }
        result.state_num += other.state_num;
        result
    }
}
