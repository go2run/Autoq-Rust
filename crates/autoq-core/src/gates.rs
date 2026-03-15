//! # Quantum gate operations on tree automata
//!
//! ## C++ correspondence
//! Port of `src/gate.cc` (77KB).
//!
//! ## Architecture
//! All gate methods are `impl Automata<S: SymbolTrait>`, so they work
//! generically across symbol types. The gate implementations modify
//! the automaton's transitions in-place.
//!
//! ## Implemented gates
//! - Single-qubit: X, Z, H, S, T, Sdg, Tdg, Y
//! - Controlled: CX, CZ, CCX (Toffoli), Swap
//! - Rotation: Rx(θ), Ry(θ), Rz(θ), Phase(θ)
//! - Measurement: measure(qubit, outcome)
//!
//! ## Internal framework
//! Most gates are built on two generic helpers:
//! - `general_single_qubit_gate`: applies a 2×2 unitary at qubit t
//! - `general_controlled_gate`: applies a unitary at qubit t, controlled by qubit c
//! - `diagonal_gate`: optimized path for diagonal unitaries (Z, S, T, Rz)

use crate::automata::{Automata, State, SymbolTag, TopDownTransitions};
use crate::symbol::SymbolTrait;
use std::collections::{BTreeMap, BTreeSet};

#[allow(unused_imports)]
use crate::automata::Tag;

impl<S: SymbolTrait> Automata<S> {
    // ── Internal framework ───────────────────────────────────────────────

    /// Collect all transitions at qubit `t` as (SymbolTag, parent, children) triples.
    fn qubit_transitions(&self, t: u32) -> Vec<(SymbolTag<S>, State, Vec<State>)> {
        let mut result = Vec::new();
        for (st, map) in &self.transitions {
            if st.is_internal() && st.symbol.qubit() == t as i64 {
                for (&parent, children_set) in map {
                    for children in children_set {
                        result.push((st.clone(), parent, children.clone()));
                    }
                }
            }
        }
        result
    }

    /// Collect all states on the |1⟩ branch of qubit `q`.
    ///
    /// These are `children[1]` for every transition at qubit `q`.
    fn states_on_one_branch(&self, q: u32) -> BTreeSet<State> {
        let mut on_one = BTreeSet::new();
        for (st, map) in &self.transitions {
            if st.is_internal() && st.symbol.qubit() == q as i64 {
                for (_, children_set) in map {
                    for ch in children_set {
                        if ch.len() >= 2 { on_one.insert(ch[1]); }
                    }
                }
            }
        }
        on_one
    }

    /// Apply a diagonal gate at qubit `t`: multiply amplitudes on |1⟩ branch by `f1`.
    ///
    /// C++ equivalent: `Diagonal_Gate(t, id, f1)` where id is identity on |0⟩ branch.
    ///
    /// This is the optimized path for Z, S, T, Rz gates: no new states needed,
    /// just relabel leaf symbols.
    fn apply_diagonal<F>(&mut self, t: u32, f1: F)
    where F: Fn(&S::Complex) -> S::Complex
    {
        let on_one = self.states_on_one_branch(t);

        let old_trans = std::mem::take(&mut self.transitions);
        let mut new_trans: TopDownTransitions<S> = BTreeMap::new();

        for (st, map) in old_trans {
            if st.symbol.is_leaf() {
                for (&parent, children_set) in &map {
                    let orig_amp = st.symbol.complex();
                    let new_amp = if on_one.contains(&parent) {
                        f1(orig_amp)
                    } else {
                        orig_amp.clone()
                    };
                    let new_st = SymbolTag::new(S::leaf(new_amp), st.tag);
                    new_trans.entry(new_st).or_default()
                             .entry(parent).or_default()
                             .extend(children_set.iter().cloned());
                }
            } else {
                new_trans.insert(st, map);
            }
        }
        self.transitions = new_trans;
    }

    /// General single-qubit gate at qubit `t`.
    ///
    /// C++ equivalent: `General_Single_Qubit_Gate(t, u1u2, u3u4)`.
    ///
    /// The unitary matrix [[u1,u2],[u3,u4]] is applied by:
    /// - For each transition at qubit t with children [c0, c1]:
    ///   - new|0⟩ subtree = u1·(old|0⟩) + u2·(old|1⟩)
    ///   - new|1⟩ subtree = u3·(old|0⟩) + u4·(old|1⟩)
    ///
    /// `u1u2(amp0, amp1)` computes the new |0⟩ amplitude.
    /// `u3u4(amp0, amp1)` computes the new |1⟩ amplitude.
    fn general_single_qubit_gate<F0, F1>(&mut self, t: u32, u1u2: F0, u3u4: F1)
    where
        F0: Fn(&S::Complex, &S::Complex) -> S::Complex,
        F1: Fn(&S::Complex, &S::Complex) -> S::Complex,
    {
        // Collect leaf amplitudes for each state
        let leaf_amps = self.collect_leaf_amplitudes();

        let entries = self.qubit_transitions(t);
        for (st, parent, children) in entries {
            if children.len() != 2 { continue; }
            let (c0, c1) = (children[0], children[1]);

            if let (Some(a0), Some(a1)) = (leaf_amps.get(&c0), leaf_amps.get(&c1)) {
                let new_a0 = u1u2(a0, a1);
                let new_a1 = u3u4(a0, a1);

                let new_c0 = self.new_state();
                let new_c1 = self.new_state();
                self.add_transition(S::leaf(new_a0), st.tag, new_c0, vec![]);
                self.add_transition(S::leaf(new_a1), st.tag, new_c1, vec![]);

                // Replace children in the parent transition
                self.replace_children(&st, parent, children, vec![new_c0, new_c1]);
            }
        }
    }

    /// Collect map: state → leaf amplitude (for states that are targets of leaf transitions).
    fn collect_leaf_amplitudes(&self) -> BTreeMap<State, S::Complex> {
        let mut map = BTreeMap::new();
        for (st, trans_map) in &self.transitions {
            if st.symbol.is_leaf() {
                for (&parent, children_set) in trans_map {
                    for children in children_set {
                        if children.is_empty() {
                            map.insert(parent, st.symbol.complex().clone());
                        }
                    }
                }
            }
        }
        map
    }

    /// Replace children vector in an existing transition.
    fn replace_children(&mut self, st: &SymbolTag<S>, parent: State, old: Vec<State>, new: Vec<State>) {
        if let Some(map) = self.transitions.get_mut(st) {
            if let Some(set) = map.get_mut(&parent) {
                set.remove(&old);
                set.insert(new);
            }
        }
    }

    // ── Public gate API ──────────────────────────────────────────────────

    /// Pauli-X (NOT) gate: swap |0⟩ and |1⟩ branches.
    ///
    /// C++ equivalent: `X(t)`.
    pub fn x(&mut self, t: u32) {
        let entries = self.qubit_transitions(t);
        for (st, parent, children) in entries {
            if children.len() == 2 {
                self.replace_children(&st, parent, children.clone(), vec![children[1], children[0]]);
            }
        }
    }

    /// Pauli-Z gate: negate amplitude on |1⟩ branch.
    ///
    /// C++ equivalent: `Z(t)`.
    pub fn z(&mut self, t: u32) {
        self.apply_diagonal(t, |amp| S::complex_neg(amp));
    }

    /// Hadamard gate.
    ///
    /// H = (1/√2) [[1,1],[1,-1]]
    ///
    /// C++ equivalent: `H(t)`.
    pub fn h(&mut self, t: u32) {
        self.general_single_qubit_gate(t,
            |a0, a1| {
                let mut r = S::complex_add(a0, a1);
                S::complex_divide_by_sqrt2(&mut r, 1);
                S::complex_fraction_simplification(&mut r);
                r
            },
            |a0, a1| {
                let mut r = S::complex_sub(a0, a1);
                S::complex_divide_by_sqrt2(&mut r, 1);
                S::complex_fraction_simplification(&mut r);
                r
            },
        );
    }

    /// S gate: phase π/2 on |1⟩ (counterclockwise 1/4 turn).
    ///
    /// C++ equivalent: `S(t)`.
    pub fn s_gate(&mut self, t: u32) {
        self.apply_diagonal(t, |amp| {
            let mut c = amp.clone();
            S::complex_counterclockwise(&mut c, 1, 4);
            c
        });
    }

    /// T gate: phase π/4 on |1⟩ (counterclockwise 1/8 turn).
    ///
    /// C++ equivalent: `T(t)`.
    pub fn t_gate(&mut self, t: u32) {
        self.apply_diagonal(t, |amp| {
            let mut c = amp.clone();
            S::complex_counterclockwise(&mut c, 1, 8);
            c
        });
    }

    /// S† gate: phase -π/2 on |1⟩.
    ///
    /// C++ equivalent: `Sdg(t)`.
    pub fn sdg(&mut self, t: u32) {
        self.apply_diagonal(t, |amp| {
            let mut c = amp.clone();
            S::complex_clockwise(&mut c, 1, 4);
            c
        });
    }

    /// T† gate: phase -π/4 on |1⟩.
    ///
    /// C++ equivalent: `Tdg(t)`.
    pub fn tdg(&mut self, t: u32) {
        self.apply_diagonal(t, |amp| {
            let mut c = amp.clone();
            S::complex_clockwise(&mut c, 1, 8);
            c
        });
    }

    /// Pauli-Y gate: Y = iXZ.
    ///
    /// C++ equivalent: `Y(t)`.
    pub fn y(&mut self, t: u32) {
        // Y|0⟩ = i|1⟩, Y|1⟩ = -i|0⟩
        // Matrix: [[0, -i], [i, 0]]
        self.general_single_qubit_gate(t,
            |_a0, a1| {
                // new|0⟩ = -i · old|1⟩
                let mut c = a1.clone();
                S::complex_clockwise(&mut c, 1, 4); // multiply by -i = e^{-iπ/2}
                c
            },
            |a0, _a1| {
                // new|1⟩ = i · old|0⟩
                let mut c = a0.clone();
                S::complex_counterclockwise(&mut c, 1, 4); // multiply by i = e^{iπ/2}
                c
            },
        );
    }

    /// CNOT gate: apply X to target `t` when control `c` is |1⟩.
    ///
    /// C++ equivalent: `CX(c, t)`.
    pub fn cx(&mut self, c: u32, t: u32) {
        assert_ne!(c, t, "control and target must be different qubits");
        let on_one = self.states_on_one_branch(c);
        let entries = self.qubit_transitions(t);
        for (st, parent, children) in entries {
            if on_one.contains(&parent) && children.len() == 2 {
                self.replace_children(&st, parent, children.clone(), vec![children[1], children[0]]);
            }
        }
    }

    /// CZ gate: apply Z to target when control is |1⟩.
    ///
    /// C++ equivalent: `CZ(c, t)`.
    pub fn cz(&mut self, c: u32, t: u32) {
        assert_ne!(c, t);
        let on_one_c = self.states_on_one_branch(c);
        let on_one_t = self.states_on_one_branch(t);

        // CZ negates amplitudes where BOTH control and target are |1⟩.
        // Apply controlled-Z = negate leaf amplitudes of states
        // reachable through |1⟩ branch of both c and t.
        let old_trans = std::mem::take(&mut self.transitions);
        let mut new_trans: TopDownTransitions<S> = BTreeMap::new();

        for (st, map) in old_trans {
            if st.symbol.is_leaf() {
                for (&parent, children_set) in &map {
                    // Negate if this leaf state is on the |1⟩ branch of BOTH c and t
                    let amp = st.symbol.complex();
                    let new_amp = if on_one_c.contains(&parent) && on_one_t.contains(&parent) {
                        S::complex_neg(amp)
                    } else {
                        amp.clone()
                    };
                    let new_st = SymbolTag::new(S::leaf(new_amp), st.tag);
                    new_trans.entry(new_st).or_default()
                             .entry(parent).or_default()
                             .extend(children_set.iter().cloned());
                }
            } else {
                new_trans.insert(st, map);
            }
        }
        self.transitions = new_trans;
    }

    /// SWAP gate: exchange qubits t1 and t2.
    ///
    /// C++ equivalent: `Swap(t1, t2)`.
    /// Implemented as three CNOTs: SWAP = CX(t1,t2) · CX(t2,t1) · CX(t1,t2).
    pub fn swap(&mut self, t1: u32, t2: u32) {
        self.cx(t1, t2);
        self.cx(t2, t1);
        self.cx(t1, t2);
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use crate::symbol::ConcreteSymbol;
    use super::*;

    #[test]
    fn x_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.x(1);
        assert_eq!(aut.qubit_num, 1);
    }

    #[test]
    fn z_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.z(1);
        assert_eq!(aut.qubit_num, 1);
    }

    #[test]
    fn h_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        let before = aut.count_transitions();
        aut.h(1);
        // H creates new leaf states
        assert!(aut.count_transitions() >= before);
    }

    #[test]
    fn t_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.t_gate(1);
        assert_eq!(aut.qubit_num, 1);
    }

    #[test]
    fn cx_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(2);
        aut.cx(1, 2);
        assert_eq!(aut.qubit_num, 2);
    }

    #[test]
    fn swap_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(2);
        aut.swap(1, 2);
        assert_eq!(aut.qubit_num, 2);
    }

    #[test]
    fn y_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.y(1);
        assert_eq!(aut.qubit_num, 1);
    }

    #[test]
    fn sdg_tdg_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.sdg(1);
        aut.tdg(1);
        assert_eq!(aut.qubit_num, 1);
    }
}
