//! # Quantum gate operations on tree automata
//!
//! ## C++ correspondence
//! Faithful port of `src/gate.cc` (77KB).
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
//! Most gates are built on three generic helpers:
//! - `general_single_qubit_gate`: product construction for 2×2 unitary at qubit t
//! - `general_controlled_gate`: product construction with control qubit(s)
//! - `diagonal_gate`: optimized tree-splitting for diagonal unitaries (Z, S, T, Rz)

use crate::automata::{Automata, State, SymbolTag, Tag, TopDownTransitions};
use crate::symbol::SymbolTrait;
use std::collections::{BTreeMap, BTreeSet};

/// Internal transitions indexed by qubit level.
/// `internal_transitions[q]` = map from Tag → (parent_state → set of child vectors)
/// C++ equivalent: `InternalTopDownTransitions` = `vector<map<Tag, map<State, set<StateVector>>>>`.
type InternalTopDownTransitions = Vec<BTreeMap<Tag, BTreeMap<State, BTreeSet<Vec<State>>>>>;

impl<S: SymbolTrait> Automata<S> {
    // ── State-pair encoding (C++ macros L/R) ────────────────────────────

    /// Encode left state pair: `stateNum + s1 * stateNum + s2`
    #[inline]
    fn l(&self, s1: State, s2: State) -> State {
        self.state_num + s1 * self.state_num + s2
    }

    /// Encode right state pair: `stateNum + stateNum^2 + s1 * stateNum + s2`
    #[inline]
    fn r(&self, s1: State, s2: State) -> State {
        self.state_num + self.state_num * self.state_num + s1 * self.state_num + s2
    }

    // ── General Single Qubit Gate (faithful port of C++) ────────────────

    /// General single-qubit gate at qubit `t` using product construction.
    ///
    /// C++ equivalent: `General_Single_Qubit_Gate(t, u1u2, u3u4)`.
    ///
    /// The unitary matrix [[u1,u2],[u3,u4]] is applied via product construction:
    /// - States are encoded as pairs L(s1,s2) and R(s1,s2)
    /// - L represents |0⟩ component, R represents |1⟩ component
    /// - Leaf amplitudes are combined via u1u2 and u3u4 functions
    fn general_single_qubit_gate<F0, F1>(&mut self, t: i64, u1u2: F0, u3u4: F1)
    where
        F0: Fn(&S, &S) -> S,
        F1: Fn(&S, &S) -> S,
    {
        let sn = self.state_num;
        let mut result = self.clone_metadata();
        let max_paired = self.r(sn - 1, sn - 1) + 1;

        // We assume transitions are ordered by symbols (BTreeMap guarantees this).
        // Phase 1: Copy transitions for qubits < t
        for (st, map) in &self.transitions {
            if st.is_internal() && st.symbol.qubit() < t {
                result.transitions.insert(st.clone(), map.clone());
            }
        }

        // Phase 2: Split transitions at qubit == t into L/R pairs
        let mut possible_next = vec![false; max_paired as usize];
        for (st, map) in &self.transitions {
            if !st.is_internal() || st.symbol.qubit() != t { continue; }
            let ref1 = result.transitions.entry(st.clone()).or_default();
            for (&parent, children_set) in map {
                let ref2 = ref1.entry(parent).or_default();
                for children in children_set {
                    assert_eq!(children.len(), 2);
                    let l = self.l(children[0], children[1]);
                    let r = self.r(children[0], children[1]);
                    ref2.insert(vec![l, r]);
                    possible_next[l as usize] = true;
                    possible_next[r as usize] = true;
                }
            }
        }

        // Phase 3: Internal transitions for qubits > t (product construction)
        // qcfi: state → tag → symbol → vec of child-state-vectors
        let mut qcfi: BTreeMap<State, BTreeMap<Tag, BTreeMap<S, Vec<Vec<State>>>>> = BTreeMap::new();
        let mut possible_prev = possible_next.clone();
        let mut prev_qubit: Option<i64> = None;
        let mut layer_entries: Vec<(&SymbolTag<S>, &BTreeMap<State, BTreeSet<Vec<State>>>)> = Vec::new();

        // Collect all internal transitions with qubit > t
        let internal_gt: Vec<_> = self.transitions.iter()
            .filter(|(st, _)| st.is_internal() && st.symbol.qubit() > t)
            .collect();

        for (st, map) in &internal_gt {
            let qubit = st.symbol.qubit();
            if prev_qubit.is_some() && prev_qubit != Some(qubit) {
                // Flush previous layer
                Self::flush_qcfi(&qcfi, &mut result.transitions, &mut possible_next);
                qcfi.clear();
                possible_prev = possible_next.clone();
            }
            if prev_qubit != Some(qubit) {
                layer_entries.clear();
                // Collect all entries at this qubit level
                for (st2, map2) in &internal_gt {
                    if st2.symbol.qubit() == qubit {
                        layer_entries.push((st2, map2));
                    }
                }
            }
            prev_qubit = Some(qubit);

            for (st2, map2) in &layer_entries {
                let color_intersection = st.tag & st2.tag;
                if color_intersection == 0 { continue; }
                for (&top1, ins1) in *map {
                    for (&top2, ins2) in *map2 {
                        let l_top = self.l(top1, top2);
                        if self.has_loop || (l_top < possible_prev.len() as State && possible_prev[l_top as usize]) {
                            let ref_l = qcfi.entry(l_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(st.symbol.clone()).or_default();
                            for in1 in ins1 {
                                for in2 in ins2 {
                                    ref_l.push(vec![self.l(in1[0], in2[0]), self.l(in1[1], in2[1])]);
                                }
                            }
                        }
                        let r_top = self.r(top1, top2);
                        if self.has_loop || (r_top < possible_prev.len() as State && possible_prev[r_top as usize]) {
                            let ref_r = qcfi.entry(r_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(st.symbol.clone()).or_default();
                            for in1 in ins1 {
                                for in2 in ins2 {
                                    ref_r.push(vec![self.r(in1[0], in2[0]), self.r(in1[1], in2[1])]);
                                }
                            }
                        }
                    }
                }
            }
        }
        // Flush last internal layer
        if !qcfi.is_empty() {
            Self::flush_qcfi(&qcfi, &mut result.transitions, &mut possible_next);
            qcfi.clear();
        }

        // Phase 4: Leaf transitions (apply unitary)
        possible_prev = possible_next;
        qcfi.clear();
        let leaf_entries: Vec<_> = self.transitions.iter()
            .filter(|(st, _)| st.is_leaf())
            .collect();

        for (st, map) in &leaf_entries {
            for (st2, map2) in &leaf_entries {
                let color_intersection = st.tag & st2.tag;
                if color_intersection == 0 { continue; }
                for (&top1, _) in *map {
                    for (&top2, _) in *map2 {
                        let l_top = self.l(top1, top2);
                        if self.has_loop || (l_top < possible_prev.len() as State && possible_prev[l_top as usize]) {
                            let new_sym = u1u2(&st.symbol, &st2.symbol);
                            qcfi.entry(l_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(new_sym).or_default()
                                .push(vec![]);
                        }
                        let r_top = self.r(top1, top2);
                        if self.has_loop || (r_top < possible_prev.len() as State && possible_prev[r_top as usize]) {
                            let new_sym = u3u4(&st.symbol, &st2.symbol);
                            qcfi.entry(r_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(new_sym).or_default()
                                .push(vec![]);
                        }
                    }
                }
            }
        }
        for (q, tag_map) in &qcfi {
            for (&tag, sym_map) in tag_map {
                for (sym, vecs) in sym_map {
                    let st = SymbolTag::new(sym.clone(), tag);
                    let entry = result.transitions.entry(st).or_default()
                        .entry(*q).or_default();
                    for v in vecs {
                        entry.insert(v.clone());
                    }
                }
            }
        }

        result.state_num = self.r(sn - 1, sn - 1) + 1;
        result.reduce();
        *self = result;
    }

    /// Flush qcfi accumulator into result transitions and update possible_next_level_states.
    fn flush_qcfi(
        qcfi: &BTreeMap<State, BTreeMap<Tag, BTreeMap<S, Vec<Vec<State>>>>>,
        result_trans: &mut TopDownTransitions<S>,
        possible_next: &mut Vec<bool>,
    ) {
        for (&q, tag_map) in qcfi {
            for (&tag, sym_map) in tag_map {
                for (sym, vecs) in sym_map {
                    let st = SymbolTag::new(sym.clone(), tag);
                    let entry = result_trans.entry(st).or_default()
                        .entry(q).or_default();
                    for v in vecs {
                        entry.insert(v.clone());
                        for &s in v {
                            if (s as usize) < possible_next.len() {
                                possible_next[s as usize] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // ── Diagonal Gate (faithful port of C++) ────────────────────────────

    /// Diagonal gate at qubit `t`: multiply leaf amplitudes by c0 (|0⟩ branch) and c1 (|1⟩ branch).
    ///
    /// C++ equivalent: `Diagonal_Gate(t, multiply_by_c0, multiply_by_c1)`.
    ///
    /// Uses tree-splitting: states below qubit t are duplicated into "left" (|0⟩) and
    /// "right" (|1⟩) copies. Each copy gets its leaf amplitudes multiplied by the
    /// corresponding function.
    fn diagonal_gate<F0, F1>(&mut self, t: i64, multiply_by_c0: F0, multiply_by_c1: F1)
    where
        F0: Fn(&mut S),
        F1: Fn(&mut S),
    {
        let sn = self.state_num;
        let mut transitions2: TopDownTransitions<S> = BTreeMap::new();
        let mut top_lr: BTreeMap<State, u8> = BTreeMap::new();
        let mut child_lr: BTreeMap<State, u8> = BTreeMap::new();

        // Convert to InternalTopDownTransitions (indexed by qubit level)
        let loop_extra = if self.has_loop { 1 } else { 0 };
        let mut internal_trans: InternalTopDownTransitions =
            vec![BTreeMap::new(); (self.qubit_num + 1 + loop_extra) as usize];
        let mut leaf_trans: TopDownTransitions<S> = BTreeMap::new();

        for (st, map) in &self.transitions {
            if st.is_internal() {
                let q = st.symbol.qubit();
                if q < t {
                    transitions2.insert(st.clone(), map.clone());
                } else {
                    internal_trans[q as usize].entry(st.tag).or_default()
                        .extend(map.iter().map(|(&k, v)| (k, v.clone())));
                }
            } else {
                for (&parent, children_set) in map {
                    leaf_trans.entry(st.clone()).or_default()
                        .entry(parent).or_default()
                        .extend(children_set.iter().cloned());
                }
            }
        }

        let max_q = self.qubit_num as i64 + loop_extra as i64;
        for q in t..=max_q {
            if q == t {
                // Construct childStateIsLeftOrRight
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (_, children_set) in out_ins {
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= 0b10; // left (|0⟩)
                                *child_lr.entry(children[1]).or_insert(0) |= 0b01; // right (|1⟩)
                            }
                        }
                    }
                }
                // Emit transitions at level t: split right child
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&parent, children_set) in out_ins {
                            let ref_parent = ref_map.entry(parent).or_default();
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                let new_in1 = children[1] + sn; // queryChildID
                                ref_parent.insert(vec![children[0], new_in1]);
                            }
                        }
                    }
                }
            } else {
                // q > t: propagate left/right info
                // First, collect parent states at this level
                let mut parents_at_level: BTreeSet<State> = BTreeSet::new();
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, _) in out_ins {
                            parents_at_level.insert(top);
                        }
                    }
                }
                // Carry forward markings for non-parent states (leaf states that
                // skip this level). Without this, states marked at level t but not
                // appearing as parents at subsequent levels lose their markings.
                for (&state, &val) in &top_lr {
                    if !parents_at_level.contains(&state) {
                        *child_lr.entry(state).or_insert(0) |= val;
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                            }
                        }
                    }
                }
                // Emit transitions
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                // Original tree
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            let new_top = top + sn; // queryTopID
                            if val & 0b01 != 0 {
                                // Copied tree
                                let ref_new = ref_map.entry(new_top).or_default();
                                for children in children_set {
                                    assert_eq!(children.len(), 2);
                                    let new_in0 = children[0] + sn; // queryChildID
                                    let new_in1 = children[1] + sn;
                                    ref_new.insert(vec![new_in0, new_in1]);
                                }
                            }
                        }
                    }
                }
            }
            top_lr = child_lr;
            child_lr = BTreeMap::new();
        }

        // Leaf transitions
        for (st, map) in &leaf_trans {
            for (&top, _children_set) in map {
                let val = *top_lr.get(&top).unwrap_or(&0);
                if val & 0b10 != 0 {
                    let mut symbol_tag = st.clone();
                    multiply_by_c0(&mut symbol_tag.symbol);
                    transitions2.entry(symbol_tag).or_default()
                        .entry(top).or_default()
                        .insert(vec![]);
                }
                if val & 0b01 != 0 {
                    let mut symbol_tag = st.clone();
                    multiply_by_c1(&mut symbol_tag.symbol);
                    let new_top = top + sn; // queryTopID
                    transitions2.entry(symbol_tag).or_default()
                        .entry(new_top).or_default()
                        .insert(vec![]);
                }
            }
        }

        self.transitions = transitions2;
        self.state_num = sn * 2;
    }

    // ── General Controlled Gate (faithful port of C++) ──────────────────

    /// General controlled gate: apply unitary at qubit `t` controlled by qubit(s).
    ///
    /// C++ equivalent: `General_Controlled_Gate(c, c2, t, u1u2, u3u4, multiply_by_c0)`.
    ///
    /// Requires all control qubits > target qubit.
    /// `multiply_by_c0` transforms leaf symbols when control is |0⟩ (identity for CX).
    fn general_controlled_gate<F0, F1, FC>(
        &mut self, c: i64, c2: i64, t: i64,
        u1u2: F0, u3u4: F1, multiply_by_c0: FC,
    )
    where
        F0: Fn(&S, &S) -> S,
        F1: Fn(&S, &S) -> S,
        FC: Fn(&S) -> S,
    {
        let min_c = c.min(c2);
        assert!(min_c > t, "All control qubits must be > target qubit");

        let sn = self.state_num;
        let mut result = self.clone_metadata();
        let max_paired = self.r(sn - 1, sn - 1) + 1;

        // Phase 0: Push default leaves (operated with multiply_by_c0) and transitions < t
        for (st, map) in &self.transitions {
            if st.is_leaf() {
                let mut new_st = st.clone();
                new_st.symbol = multiply_by_c0(&st.symbol);
                for (&parent, children_set) in map {
                    result.transitions.entry(new_st.clone()).or_default()
                        .entry(parent).or_default()
                        .extend(children_set.iter().cloned());
                }
            } else if st.symbol.qubit() < t {
                result.transitions.insert(st.clone(), map.clone());
            }
        }

        // Find the first transition at qubit >= t
        let internal_ge_t: Vec<_> = self.transitions.iter()
            .filter(|(st, _)| st.is_internal() && st.symbol.qubit() >= t)
            .collect();

        // Phase 1: Split at qubit == t
        let mut possible_next = vec![false; max_paired as usize];
        for (st, map) in &internal_ge_t {
            if st.symbol.qubit() != t { continue; }
            let ref1 = result.transitions.entry((*st).clone()).or_default();
            for (&parent, children_set) in *map {
                let ref2 = ref1.entry(parent).or_default();
                for children in children_set {
                    assert_eq!(children.len(), 2);
                    let l = self.l(children[0], children[1]);
                    let r = self.r(children[0], children[1]);
                    ref2.insert(vec![l, r]);
                    possible_next[l as usize] = true;
                    possible_next[r as usize] = true;
                }
            }
        }

        // Phase 2: Internal transitions for qubits > t
        let mut qcfi: BTreeMap<State, BTreeMap<Tag, BTreeMap<S, Vec<Vec<State>>>>> = BTreeMap::new();
        let mut possible_prev = possible_next.clone();
        let mut prev_qubit: Option<i64> = None;

        let internal_gt: Vec<_> = internal_ge_t.iter()
            .filter(|(st, _)| st.symbol.qubit() > t)
            .cloned()
            .collect();

        // Collect entries by layer
        let mut layer_starts: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
        for (i, (st, _)) in internal_gt.iter().enumerate() {
            layer_starts.entry(st.symbol.qubit()).or_default().push(i);
        }

        for (st, map) in &internal_gt {
            let qubit = st.symbol.qubit();
            if prev_qubit.is_some() && prev_qubit != Some(qubit) {
                Self::flush_qcfi(&qcfi, &mut result.transitions, &mut possible_next);
                qcfi.clear();
                possible_prev = possible_next.clone();
            }
            prev_qubit = Some(qubit);

            // Get all entries at this qubit level
            let layer_indices = layer_starts.get(&qubit).unwrap();

            for &idx2 in layer_indices {
                let (st2, map2) = &internal_gt[idx2];
                let color_intersection = st.tag & st2.tag;
                if color_intersection == 0 { continue; }
                let qubit_is_not_control = qubit != c && qubit != c2;

                for (&top1, ins1) in *map {
                    for (&top2, ins2) in *map2 {
                        let l_top = self.l(top1, top2);
                        if self.has_loop || (l_top < possible_prev.len() as State && possible_prev[l_top as usize]) {
                            let ref_l = qcfi.entry(l_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(st.symbol.clone()).or_default();
                            if qubit_is_not_control {
                                for in1 in ins1 {
                                    for in2 in ins2 {
                                        ref_l.push(vec![self.l(in1[0], in2[0]), self.l(in1[1], in2[1])]);
                                    }
                                }
                            } else {
                                // Control qubit: L keeps left child from it, right child from product
                                for in1 in ins1 {
                                    for in2 in ins2 {
                                        ref_l.push(vec![in1[0], self.l(in1[1], in2[1])]);
                                    }
                                }
                            }
                        }
                        let r_top = self.r(top1, top2);
                        if self.has_loop || (r_top < possible_prev.len() as State && possible_prev[r_top as usize]) {
                            let ref_r = qcfi.entry(r_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(st.symbol.clone()).or_default();
                            if qubit_is_not_control {
                                for in1 in ins1 {
                                    for in2 in ins2 {
                                        ref_r.push(vec![self.r(in1[0], in2[0]), self.r(in1[1], in2[1])]);
                                    }
                                }
                            } else {
                                // Control qubit: R keeps right child from it2, right child from product
                                for in1 in ins1 {
                                    for in2 in ins2 {
                                        ref_r.push(vec![in2[0], self.r(in1[1], in2[1])]);
                                    }
                                }
                            }
                        }
                    }
                }

                // For qubits above min_c: also keep original transitions
                if qubit > min_c {
                    for (&top, children_set) in *map {
                        if self.has_loop || (top < possible_prev.len() as State && possible_prev[top as usize]) {
                            let ref_q = qcfi.entry(top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(st.symbol.clone()).or_default();
                            for ch in children_set {
                                ref_q.push(ch.clone());
                            }
                        }
                    }
                    for (&top, children_set) in *map2 {
                        if self.has_loop || (top < possible_prev.len() as State && possible_prev[top as usize]) {
                            let ref_q = qcfi.entry(top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(st2.symbol.clone()).or_default();
                            for ch in children_set {
                                ref_q.push(ch.clone());
                            }
                        }
                    }
                }
            }
        }
        if !qcfi.is_empty() {
            Self::flush_qcfi(&qcfi, &mut result.transitions, &mut possible_next);
            qcfi.clear();
        }

        // Phase 3: Leaf transitions (apply unitary)
        possible_prev = possible_next;
        let leaf_entries: Vec<_> = self.transitions.iter()
            .filter(|(st, _)| st.is_leaf())
            .collect();

        for (st, map) in &leaf_entries {
            for (st2, map2) in &leaf_entries {
                let color_intersection = st.tag & st2.tag;
                if color_intersection == 0 { continue; }
                for (&top1, _) in *map {
                    for (&top2, _) in *map2 {
                        let l_top = self.l(top1, top2);
                        if self.has_loop || (l_top < possible_prev.len() as State && possible_prev[l_top as usize]) {
                            let new_sym = u1u2(&st.symbol, &st2.symbol);
                            qcfi.entry(l_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(new_sym).or_default()
                                .push(vec![]);
                        }
                        let r_top = self.r(top1, top2);
                        if self.has_loop || (r_top < possible_prev.len() as State && possible_prev[r_top as usize]) {
                            let new_sym = u3u4(&st.symbol, &st2.symbol);
                            qcfi.entry(r_top).or_default()
                                .entry(color_intersection).or_default()
                                .entry(new_sym).or_default()
                                .push(vec![]);
                        }
                    }
                }
            }
        }
        for (q, tag_map) in &qcfi {
            for (&tag, sym_map) in tag_map {
                for (sym, vecs) in sym_map {
                    let st = SymbolTag::new(sym.clone(), tag);
                    let entry = result.transitions.entry(st).or_default()
                        .entry(*q).or_default();
                    for v in vecs {
                        entry.insert(v.clone());
                    }
                }
            }
        }

        result.state_num = self.r(sn - 1, sn - 1) + 1;
        *self = result;
    }

    /// Clone metadata without transitions.
    fn clone_metadata(&self) -> Self {
        Automata {
            name: self.name.clone(),
            final_states: self.final_states.clone(),
            state_num: self.state_num,
            qubit_num: self.qubit_num,
            symbolic_vars_num: self.symbolic_vars_num,
            transitions: BTreeMap::new(),
            vars: self.vars.clone(),
            constraints: self.constraints.clone(),
            has_loop: self.has_loop,
            is_topdown_deterministic: self.is_topdown_deterministic,
        }
    }

    // ── CX helper for c < t (diagonal-style splitting) ─────────────────

    /// CX implementation when control qubit < target qubit.
    ///
    /// C++ equivalent: CX(c, t) when c < t.
    ///
    /// Uses diagonal-gate-like tree splitting from the control qubit downward.
    /// At the target qubit, the |1⟩-controlled copy swaps its children.
    fn cx_c_less_than_t(&mut self, c: i64, t: i64) {
        let sn = self.state_num;
        let mut transitions2: TopDownTransitions<S> = BTreeMap::new();
        let mut top_lr: BTreeMap<State, u8> = BTreeMap::new();
        let mut child_lr: BTreeMap<State, u8> = BTreeMap::new();

        // Convert to per-qubit indexed transitions
        let mut internal_trans: InternalTopDownTransitions =
            vec![BTreeMap::new(); (self.qubit_num + 1) as usize];
        let mut leaf_trans: TopDownTransitions<S> = BTreeMap::new();

        for (st, map) in &self.transitions {
            if st.is_internal() {
                let q = st.symbol.qubit();
                if q < c {
                    transitions2.insert(st.clone(), map.clone());
                } else {
                    internal_trans[q as usize].entry(st.tag).or_default()
                        .extend(map.iter().map(|(&k, v)| (k, v.clone())));
                }
            } else {
                for (&parent, children_set) in map {
                    leaf_trans.entry(st.clone()).or_default()
                        .entry(parent).or_default()
                        .extend(children_set.iter().cloned());
                }
            }
        }

        for q in c..=self.qubit_num as i64 {
            if q == c {
                // At control qubit: split left (|0⟩) and right (|1⟩) branches
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (_, children_set) in out_ins {
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= 0b10;
                                *child_lr.entry(children[1]).or_insert(0) |= 0b01;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&parent, children_set) in out_ins {
                            let ref_parent = ref_map.entry(parent).or_default();
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                let new_in1 = children[1] + sn;
                                ref_parent.insert(vec![children[0], new_in1]);
                            }
                        }
                    }
                }
            } else {
                // q > c: propagate left/right info
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                            }
                        }
                    }
                }
                // Emit transitions
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                // Original tree: keep as-is
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                // Copied tree: at target qubit, swap children
                                let new_top = top + sn;
                                let ref_new = ref_map.entry(new_top).or_default();
                                for children in children_set {
                                    assert_eq!(children.len(), 2);
                                    let new_in0 = children[0] + sn;
                                    let new_in1 = children[1] + sn;
                                    if q == t {
                                        // Apply X: swap children
                                        ref_new.insert(vec![new_in1, new_in0]);
                                    } else {
                                        ref_new.insert(vec![new_in0, new_in1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            top_lr = child_lr;
            child_lr = BTreeMap::new();
        }

        // Leaf transitions
        for (st, map) in &leaf_trans {
            let ref_map = transitions2.entry(st.clone()).or_default();
            for (&top, _) in map {
                let val = *top_lr.get(&top).unwrap_or(&0);
                if val & 0b10 != 0 {
                    ref_map.entry(top).or_default().insert(vec![]);
                }
                if val & 0b01 != 0 {
                    let new_top = top + sn;
                    ref_map.entry(new_top).or_default().insert(vec![]);
                }
            }
        }

        self.transitions = transitions2;
        self.state_num = sn * 2;
    }

    // ── CZ (faithful port of C++) ──────────────────────────────────────

    /// CZ implementation: triple state space splitting.
    ///
    /// C++ equivalent: `CZ(c, t)`.
    /// Ensures c < t (swaps if needed since CZ is symmetric).
    fn cz_full(&mut self, c_in: i64, t_in: i64) {
        let (c, t) = if c_in < t_in { (c_in, t_in) } else { (t_in, c_in) };
        let sn = self.state_num;
        let mut transitions2: TopDownTransitions<S> = BTreeMap::new();
        let mut top_lr: BTreeMap<State, u8> = BTreeMap::new();
        let mut child_lr: BTreeMap<State, u8> = BTreeMap::new();
        let mut top_lr2: BTreeMap<State, u8> = BTreeMap::new();
        let mut child_lr2: BTreeMap<State, u8> = BTreeMap::new();

        let mut internal_trans: InternalTopDownTransitions =
            vec![BTreeMap::new(); (self.qubit_num + 1) as usize];
        let mut leaf_trans: TopDownTransitions<S> = BTreeMap::new();

        for (st, map) in &self.transitions {
            if st.is_internal() {
                let q = st.symbol.qubit();
                if q < c {
                    transitions2.insert(st.clone(), map.clone());
                } else {
                    internal_trans[q as usize].entry(st.tag).or_default()
                        .extend(map.iter().map(|(&k, v)| (k, v.clone())));
                }
            } else {
                for (&parent, children_set) in map {
                    leaf_trans.entry(st.clone()).or_default()
                        .entry(parent).or_default()
                        .extend(children_set.iter().cloned());
                }
            }
        }

        for q in c..=self.qubit_num as i64 {
            if q == c {
                // Split at first control
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (_, children_set) in out_ins {
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= 0b10;
                                *child_lr.entry(children[1]).or_insert(0) |= 0b01;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&parent, children_set) in out_ins {
                            let ref_parent = ref_map.entry(parent).or_default();
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                let new_in1 = children[1] + sn;
                                ref_parent.insert(vec![children[0], new_in1]);
                            }
                        }
                    }
                }
            } else if q < t {
                // Between c and t: propagate first split
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                let new_top = top + sn;
                                let ref_new = ref_map.entry(new_top).or_default();
                                for children in children_set {
                                    assert_eq!(children.len(), 2);
                                    ref_new.insert(vec![children[0] + sn, children[1] + sn]);
                                }
                            }
                        }
                    }
                }
            } else if q == t {
                // At target qubit: second split for |1⟩-controlled copy
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                                if val & 0b01 != 0 {
                                    *child_lr2.entry(children[0]).or_insert(0) |= 0b10;
                                    *child_lr2.entry(children[1]).or_insert(0) |= 0b01;
                                }
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                let new_top = top + sn;
                                let ref_new = ref_map.entry(new_top).or_default();
                                for children in children_set {
                                    assert_eq!(children.len(), 2);
                                    let new_in0 = children[0] + sn;
                                    let new_in1 = children[1] + sn * 2; // queryChildID2
                                    ref_new.insert(vec![new_in0, new_in1]);
                                }
                            }
                        }
                    }
                }
            } else {
                // q > t: propagate both splits
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            let val2 = *top_lr2.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                                *child_lr2.entry(children[0]).or_insert(0) |= val2;
                                *child_lr2.entry(children[1]).or_insert(0) |= val2;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                let val2 = *top_lr2.get(&top).unwrap_or(&0);
                                if val2 & 0b10 != 0 {
                                    let new_top = top + sn;
                                    let ref_new = ref_map.entry(new_top).or_default();
                                    for children in children_set {
                                        assert_eq!(children.len(), 2);
                                        ref_new.insert(vec![children[0] + sn, children[1] + sn]);
                                    }
                                }
                                if val2 & 0b01 != 0 {
                                    let new_top = top + sn * 2; // queryTopID2
                                    let ref_new = ref_map.entry(new_top).or_default();
                                    for children in children_set {
                                        assert_eq!(children.len(), 2);
                                        ref_new.insert(vec![children[0] + sn * 2, children[1] + sn * 2]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            top_lr = child_lr;
            child_lr = BTreeMap::new();
            top_lr2 = child_lr2;
            child_lr2 = BTreeMap::new();
        }

        // Leaf transitions - collect all insertions first, then apply
        let mut leaf_insertions: Vec<(SymbolTag<S>, State)> = Vec::new();
        for (st, map) in &leaf_trans {
            for (&top, _) in map {
                let val = *top_lr.get(&top).unwrap_or(&0);
                if val & 0b10 != 0 {
                    leaf_insertions.push((st.clone(), top));
                }
                if val & 0b01 != 0 {
                    let val2 = *top_lr2.get(&top).unwrap_or(&0);
                    if val2 & 0b10 != 0 {
                        leaf_insertions.push((st.clone(), top + sn));
                    }
                    if val2 & 0b01 != 0 {
                        // Negate the amplitude for CZ
                        let negated_st = SymbolTag::new(
                            S::leaf(S::complex_neg(st.symbol.complex())),
                            st.tag,
                        );
                        leaf_insertions.push((negated_st, top + sn * 2));
                    }
                }
            }
        }
        for (st, parent) in leaf_insertions {
            transitions2.entry(st).or_default()
                .entry(parent).or_default()
                .insert(vec![]);
        }

        self.transitions = transitions2;
        self.state_num = sn * 3;
        self.reduce();
    }

    // ── CCX (faithful port of C++) ─────────────────────────────────────

    /// CCX when c < c2 < t: triple state space splitting (like CZ but with X at target).
    fn ccx_c_c2_t(&mut self, c: i64, c2: i64, t: i64) {
        let sn = self.state_num;
        let mut transitions2: TopDownTransitions<S> = BTreeMap::new();
        let mut top_lr: BTreeMap<State, u8> = BTreeMap::new();
        let mut child_lr: BTreeMap<State, u8> = BTreeMap::new();
        let mut top_lr2: BTreeMap<State, u8> = BTreeMap::new();
        let mut child_lr2: BTreeMap<State, u8> = BTreeMap::new();

        let mut internal_trans: InternalTopDownTransitions =
            vec![BTreeMap::new(); (self.qubit_num + 1) as usize];
        let mut leaf_trans: TopDownTransitions<S> = BTreeMap::new();

        for (st, map) in &self.transitions {
            if st.is_internal() {
                let q = st.symbol.qubit();
                if q < c {
                    transitions2.insert(st.clone(), map.clone());
                } else {
                    internal_trans[q as usize].entry(st.tag).or_default()
                        .extend(map.iter().map(|(&k, v)| (k, v.clone())));
                }
            } else {
                for (&parent, children_set) in map {
                    leaf_trans.entry(st.clone()).or_default()
                        .entry(parent).or_default()
                        .extend(children_set.iter().cloned());
                }
            }
        }

        for q in c..=self.qubit_num as i64 {
            if q == c {
                // First control: split left/right
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (_, children_set) in out_ins {
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= 0b10;
                                *child_lr.entry(children[1]).or_insert(0) |= 0b01;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&parent, children_set) in out_ins {
                            let ref_parent = ref_map.entry(parent).or_default();
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                ref_parent.insert(vec![children[0], children[1] + sn]);
                            }
                        }
                    }
                }
            } else if q < c2 {
                // Between c and c2: propagate first split
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                let new_top = top + sn;
                                let ref_new = ref_map.entry(new_top).or_default();
                                for children in children_set {
                                    assert_eq!(children.len(), 2);
                                    ref_new.insert(vec![children[0] + sn, children[1] + sn]);
                                }
                            }
                        }
                    }
                }
            } else if q == c2 {
                // Second control: second split
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                                if val & 0b01 != 0 {
                                    *child_lr2.entry(children[0]).or_insert(0) |= 0b10;
                                    *child_lr2.entry(children[1]).or_insert(0) |= 0b01;
                                }
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                let new_top = top + sn;
                                let ref_new = ref_map.entry(new_top).or_default();
                                for children in children_set {
                                    assert_eq!(children.len(), 2);
                                    ref_new.insert(vec![children[0] + sn, children[1] + sn * 2]);
                                }
                            }
                        }
                    }
                }
            } else {
                // q > c2: propagate both splits
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (_, out_ins) in tag_outins {
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            let val2 = *top_lr2.get(&top).unwrap_or(&0);
                            for children in children_set {
                                assert_eq!(children.len(), 2);
                                *child_lr.entry(children[0]).or_insert(0) |= val;
                                *child_lr.entry(children[1]).or_insert(0) |= val;
                                *child_lr2.entry(children[0]).or_insert(0) |= val2;
                                *child_lr2.entry(children[1]).or_insert(0) |= val2;
                            }
                        }
                    }
                }
                if let Some(tag_outins) = internal_trans.get(q as usize) {
                    for (&tag, out_ins) in tag_outins {
                        let st = SymbolTag::new(S::internal(q), tag);
                        let ref_map = transitions2.entry(st).or_default();
                        for (&top, children_set) in out_ins {
                            let val = *top_lr.get(&top).unwrap_or(&0);
                            if val & 0b10 != 0 {
                                ref_map.entry(top).or_default()
                                    .extend(children_set.iter().cloned());
                            }
                            if val & 0b01 != 0 {
                                let val2 = *top_lr2.get(&top).unwrap_or(&0);
                                if val2 & 0b10 != 0 {
                                    let new_top = top + sn;
                                    let ref_new = ref_map.entry(new_top).or_default();
                                    for children in children_set {
                                        assert_eq!(children.len(), 2);
                                        ref_new.insert(vec![children[0] + sn, children[1] + sn]);
                                    }
                                }
                                if val2 & 0b01 != 0 {
                                    let new_top = top + sn * 2;
                                    let ref_new = ref_map.entry(new_top).or_default();
                                    for children in children_set {
                                        assert_eq!(children.len(), 2);
                                        if q == t {
                                            // Apply X: swap children in the doubly-controlled copy
                                            ref_new.insert(vec![children[1] + sn * 2, children[0] + sn * 2]);
                                        } else {
                                            ref_new.insert(vec![children[0] + sn * 2, children[1] + sn * 2]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            top_lr = child_lr;
            child_lr = BTreeMap::new();
            top_lr2 = child_lr2;
            child_lr2 = BTreeMap::new();
        }

        // Leaf transitions
        for (st, map) in &leaf_trans {
            let ref_map = transitions2.entry(st.clone()).or_default();
            for (&top, _) in map {
                let val = *top_lr.get(&top).unwrap_or(&0);
                if val & 0b10 != 0 {
                    ref_map.entry(top).or_default().insert(vec![]);
                }
                if val & 0b01 != 0 {
                    let val2 = *top_lr2.get(&top).unwrap_or(&0);
                    if val2 & 0b10 != 0 {
                        let new_top = top + sn;
                        ref_map.entry(new_top).or_default().insert(vec![]);
                    }
                    if val2 & 0b01 != 0 {
                        let new_top = top + sn * 2;
                        ref_map.entry(new_top).or_default().insert(vec![]);
                    }
                }
            }
        }

        self.transitions = transitions2;
        self.state_num = sn * 3;
    }

    /// CCX when c < t < c2: uses CX(c2, t) on a copy merged at control c.
    fn ccx_c_t_c2(&mut self, c: i64, c2: i64, t: i64) {
        let sn = self.state_num;
        let mut aut2 = self.clone();
        aut2.cx(c2 as u32, t as u32);

        // Merge aut2's transitions (shifted by sn) for qubits > c and leaves
        for (st, map) in &aut2.transitions {
            if st.is_internal() && st.symbol.qubit() <= c { continue; }
            let ref_st = self.transitions.entry(st.clone()).or_default();
            for (&parent, children_set) in map {
                let new_parent = parent + sn;
                let ref_parent = ref_st.entry(new_parent).or_default();
                for children in children_set {
                    let new_children: Vec<State> = children.iter().map(|&s| s + sn).collect();
                    ref_parent.insert(new_children);
                }
            }
        }

        // At control qubit c: redirect |1⟩ children to the shifted copy
        let qubit_c_entries: Vec<_> = self.transitions.iter()
            .filter(|(st, _)| st.is_internal() && st.symbol.qubit() == c)
            .map(|(st, map)| (st.clone(), map.clone()))
            .collect();

        for (st, _map) in qubit_c_entries {
            let ref_map = self.transitions.get_mut(&st).unwrap();
            for (_, children_set) in ref_map.iter_mut() {
                let old_set: Vec<_> = children_set.iter().cloned().collect();
                children_set.clear();
                for mut ch in old_set {
                    assert_eq!(ch.len(), 2);
                    if ch[0] < sn && ch[1] < sn {
                        ch[1] += sn; // redirect |1⟩ to shifted copy
                    }
                    children_set.insert(ch);
                }
            }
        }

        self.state_num = sn + aut2.state_num;
        self.remove_useless();
    }

    // ── Public gate API ──────────────────────────────────────────────────

    /// Pauli-X (NOT) gate: swap |0⟩ and |1⟩ branches at qubit t.
    ///
    /// C++ equivalent: `X(t)`.
    pub fn x(&mut self, t: u32) {
        let t_i64 = t as i64;
        let keys: Vec<_> = self.transitions.keys()
            .filter(|st| st.is_internal() && st.symbol.qubit() == t_i64)
            .cloned()
            .collect();

        for st in keys {
            let ref_map = self.transitions.get_mut(&st).unwrap();
            for (_, children_set) in ref_map.iter_mut() {
                let old: Vec<_> = children_set.iter().cloned().collect();
                children_set.clear();
                for ch in old {
                    assert_eq!(ch.len(), 2);
                    children_set.insert(vec![ch[1], ch[0]]);
                }
            }
        }
    }

    /// Pauli-Y gate: Y = iXZ.
    ///
    /// C++ equivalent: `Y(t)`.
    pub fn y(&mut self, t: u32) {
        // C++: X(t); Diagonal_Gate(t, degree90cw, omega_multiplication(2));
        // degree90cw = clockwise(1/4 turn) = multiply by e^{-iπ/2} = -i
        // omega_multiplication(2) = counterclockwise(2/8 turn) = multiply by e^{iπ/2} = i
        self.x(t);
        self.diagonal_gate(t as i64,
            |sym: &mut S| { S::complex_clockwise(sym.complex_mut(), 1, 4); },
            |sym: &mut S| { S::complex_counterclockwise(sym.complex_mut(), 2, 8); },
        );
        self.reduce();
    }

    /// Pauli-Z gate: negate amplitude on |1⟩ branch.
    ///
    /// C++ equivalent: `Z(t)`.
    pub fn z(&mut self, t: u32) {
        self.diagonal_gate(t as i64,
            |_: &mut S| {},
            |sym: &mut S| {
                let neg = S::complex_neg(sym.complex());
                *sym = S::leaf(neg);
            },
        );
        self.reduce();
    }

    /// Hadamard gate: H = (1/√2)[[1,1],[1,-1]].
    ///
    /// C++ equivalent: `H(t)`.
    pub fn h(&mut self, t: u32) {
        self.general_single_qubit_gate(t as i64,
            |l, r| {
                let sum = S::leaf(S::complex_add(l.complex(), r.complex()));
                let mut c = sum.into_complex();
                S::complex_divide_by_sqrt2(&mut c, 1);
                S::leaf(c)
            },
            |l, r| {
                let diff = S::leaf(S::complex_sub(l.complex(), r.complex()));
                let mut c = diff.into_complex();
                S::complex_divide_by_sqrt2(&mut c, 1);
                S::leaf(c)
            },
        );
    }

    /// S gate: phase π/2 on |1⟩.
    ///
    /// C++ equivalent: `S(t)`.
    pub fn s_gate(&mut self, t: u32) {
        self.diagonal_gate(t as i64,
            |_: &mut S| {},
            |sym: &mut S| { S::complex_counterclockwise(sym.complex_mut(), 1, 4); },
        );
        self.reduce();
    }

    /// T gate: phase π/4 on |1⟩.
    ///
    /// C++ equivalent: `T(t)`.
    pub fn t_gate(&mut self, t: u32) {
        self.diagonal_gate(t as i64,
            |_: &mut S| {},
            |sym: &mut S| { S::complex_counterclockwise(sym.complex_mut(), 1, 8); },
        );
        self.reduce();
    }

    /// S† gate: phase -π/2 on |1⟩.
    ///
    /// C++ equivalent: `Sdg(t)`.
    pub fn sdg(&mut self, t: u32) {
        self.diagonal_gate(t as i64,
            |_: &mut S| {},
            |sym: &mut S| { S::complex_clockwise(sym.complex_mut(), 1, 4); },
        );
        self.reduce();
    }

    /// T† gate: phase -π/4 on |1⟩.
    ///
    /// C++ equivalent: `Tdg(t)`.
    pub fn tdg(&mut self, t: u32) {
        self.diagonal_gate(t as i64,
            |_: &mut S| {},
            |sym: &mut S| { S::complex_clockwise(sym.complex_mut(), 1, 8); },
        );
        self.reduce();
    }

    /// Rx(θ) gate: rotation about X axis by angle π·θ_num/θ_den.
    ///
    /// C++ equivalent: `Rx(theta, t)`.
    pub fn rx(&mut self, theta_num: i64, theta_den: i64, t: u32) {
        // Rx(θ) = [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
        let tn = theta_num;
        let td = theta_den * 2;
        self.general_single_qubit_gate(t as i64,
            move |l, r| {
                let cos_part = S::complex_multiply_cos(l.complex(), tn, td);
                let isin_part = S::complex_multiply_isin(r.complex(), tn, td);
                S::leaf(S::complex_sub(&cos_part, &isin_part))
            },
            move |l, r| {
                let cos_part = S::complex_multiply_cos(r.complex(), tn, td);
                let isin_part = S::complex_multiply_isin(l.complex(), tn, td);
                S::leaf(S::complex_sub(&cos_part, &isin_part))
            },
        );
    }

    /// Ry(π/2) gate.
    ///
    /// C++ equivalent: `Ry(t)`.
    pub fn ry(&mut self, t: u32) {
        self.general_single_qubit_gate(t as i64,
            |l, r| {
                let diff = S::complex_sub(l.complex(), r.complex());
                let mut c = diff;
                S::complex_divide_by_sqrt2(&mut c, 1);
                S::leaf(c)
            },
            |l, r| {
                let sum = S::complex_add(l.complex(), r.complex());
                let mut c = sum;
                S::complex_divide_by_sqrt2(&mut c, 1);
                S::leaf(c)
            },
        );
    }

    /// Rz(θ) gate: rotation about Z axis.
    ///
    /// C++ equivalent: `Rz(theta, t)`.
    pub fn rz(&mut self, theta_num: i64, theta_den: i64, t: u32) {
        // Rz(θ) = [[e^{-iθ/2}, 0], [0, e^{iθ/2}]]
        let tn = theta_num;
        let td = theta_den * 2;
        self.diagonal_gate(t as i64,
            move |sym: &mut S| { S::complex_counterclockwise(sym.complex_mut(), -tn, td); },
            move |sym: &mut S| { S::complex_counterclockwise(sym.complex_mut(), tn, td); },
        );
        self.reduce();
    }

    /// Phase gate: multiply all leaf amplitudes by e^{2πi·r}.
    ///
    /// C++ equivalent: `Phase(r)`.
    pub fn phase(&mut self, theta_num: i64, theta_den: i64) {
        let old_trans = std::mem::take(&mut self.transitions);
        for (st, map) in old_trans {
            if st.is_internal() {
                self.transitions.insert(st, map);
            } else {
                let mut c = st.symbol.complex().clone();
                S::complex_counterclockwise(&mut c, theta_num, theta_den);
                let new_st = SymbolTag::new(S::leaf(c), st.tag);
                self.transitions.entry(new_st).or_default().extend(map);
            }
        }
    }

    /// CNOT gate: apply X to target `t` when control `c` is |1⟩.
    ///
    /// C++ equivalent: `CX(c, t)`.
    /// Handles both c > t (via General_Controlled_Gate) and c < t (via tree splitting).
    pub fn cx(&mut self, c: u32, t: u32) {
        assert_ne!(c, t, "control and target must be different qubits");
        let c_i64 = c as i64;
        let t_i64 = t as i64;
        if c_i64 > t_i64 {
            self.general_controlled_gate(c_i64, c_i64, t_i64,
                |_l, r| r.clone(),
                |l, _r| l.clone(),
                |l| l.clone(),
            );
        } else {
            self.cx_c_less_than_t(c_i64, t_i64);
        }
        self.reduce();
    }

    /// CZ gate: apply Z to target when control is |1⟩.
    ///
    /// C++ equivalent: `CZ(c, t)`.
    pub fn cz(&mut self, c: u32, t: u32) {
        assert_ne!(c, t);
        self.cz_full(c as i64, t as i64);
    }

    /// CCX (Toffoli) gate: apply X to target controlled by two qubits.
    ///
    /// C++ equivalent: `CCX(c, c2, t)`.
    pub fn ccx(&mut self, c: u32, c2: u32, t: u32) {
        assert!(c != c2 && c2 != t && t != c);
        let (mut c_i, mut c2_i) = (c as i64, c2 as i64);
        let t_i = t as i64;
        if c_i > c2_i { std::mem::swap(&mut c_i, &mut c2_i); } // ensure c < c2

        if c2_i < t_i {
            // c < c2 < t
            self.ccx_c_c2_t(c_i, c2_i, t_i);
        } else if t_i < c_i {
            // t < c < c2
            self.general_controlled_gate(c_i, c2_i, t_i,
                |_l, r| r.clone(),
                |l, _r| l.clone(),
                |l| l.clone(),
            );
        } else {
            // c < t < c2
            self.ccx_c_t_c2(c_i, c2_i, t_i);
        }
        self.reduce();
    }

    /// SWAP gate: exchange qubits t1 and t2.
    ///
    /// C++ equivalent: `Swap(t1, t2)`.
    pub fn swap(&mut self, t1: u32, t2: u32) {
        self.cx(t1, t2);
        self.cx(t2, t1);
        self.cx(t1, t2);
    }

    /// Measurement: project onto |0⟩ or |1⟩ outcome at qubit t.
    ///
    /// C++ equivalent: `measure(t, outcome)`.
    pub fn measure(&self, t: u32, outcome: bool) -> Self {
        let mut aut = self.clone();
        if outcome {
            // Keep |1⟩, zero out |0⟩
            aut.diagonal_gate(t as i64,
                |sym: &mut S| { S::complex_back_to_zero(sym.complex_mut()); },
                |_: &mut S| {},
            );
        } else {
            // Keep |0⟩, zero out |1⟩
            aut.diagonal_gate(t as i64,
                |_: &mut S| {},
                |sym: &mut S| { S::complex_back_to_zero(sym.complex_mut()); },
            );
        }
        aut.reduce();
        aut
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
    fn x_gate_twice_identity() {
        for n in [1, 2, 3] {
            let original = Automata::<ConcreteSymbol>::zero_state(n);
            let mut aut = original.clone();
            aut.x(1);
            aut.x(1);
            // After X twice, should be back to original
            assert_eq!(aut.transitions, original.transitions,
                "X² should be identity for zero_state({})", n);
        }
    }

    #[test]
    fn z_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.z(1);
        assert_eq!(aut.qubit_num, 1);
    }

    #[test]
    fn z_gate_twice_identity_debug() {
        use crate::inclusion::are_equal;
        let original = Automata::<ConcreteSymbol>::zero_state(2);
        let mut after = original.clone();
        after.z(1);
        after.print_aut("after Z(1)");
        after.z(1);
        after.print_aut("after Z(1)²");
        original.print_aut("original");
        // Check language-level equality
        let fwd = crate::inclusion::is_included_in(&original, &after);
        let bwd = crate::inclusion::is_included_in(&after, &original);
        eprintln!("original <= after: {}, after <= original: {}", fwd, bwd);
        assert!(are_equal(&original, &after),
            "Z² should be identity on zero_state(2)");
    }

    #[test]
    fn h_gate_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.h(1);
        // H|0⟩ = |+⟩ = (|0⟩ + |1⟩)/√2, should have at least 1 transition
        assert!(aut.count_transitions() > 0);
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
    fn cx_smoke_reverse() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(2);
        aut.cx(2, 1);
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

    #[test]
    fn cz_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(2);
        aut.cz(1, 2);
        assert_eq!(aut.qubit_num, 2);
    }

    #[test]
    fn ccx_smoke() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(3);
        aut.ccx(1, 2, 3);
        assert_eq!(aut.qubit_num, 3);
    }

    #[test]
    fn measure_smoke() {
        let aut = Automata::<ConcreteSymbol>::zero_state(1);
        let m0 = aut.measure(1, false);
        assert_eq!(m0.qubit_num, 1);
        let m1 = aut.measure(1, true);
        assert_eq!(m1.qubit_num, 1);
    }
}
