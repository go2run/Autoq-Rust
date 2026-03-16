//! # Automata<S> — Generic non-deterministic finite tree automaton
//!
//! ## C++ correspondence
//! Port of `include/autoq/aut_description.hh` and parts of `src/query.cc`,
//! `src/instance.cc`, `src/general.cc`.
//!
//! ## Key differences from C++
//! - `Automata<TT>` → `Automata<S: SymbolTrait>` (trait-bounded generic).
//! - `SymbolTag` is a struct with named fields instead of inheriting `std::pair`.
//! - Transition maps use `BTreeMap` for deterministic ordering (matches C++'s `std::map`).
//! - No global mutable statics (C++ uses `inline static`); statistics are passed explicitly.
//!
//! ## Transition structure
//! ```text
//! TopDownTransitions = BTreeMap<SymbolTag<S>, BTreeMap<State, BTreeSet<StateVector>>>
//! ```
//! A transition `(symbol, tag)[q1, ..., qn] → q` means:
//! - `symbol`: the label (qubit index for internal, amplitude for leaf)
//! - `tag`: color bitset — paths must share at least one common color bit
//! - `[q1, ..., qn]`: child states (empty for leaves)
//! - `q`: parent state

use crate::symbol::SymbolTrait;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// State identifier (matches C++ `int64_t`).
pub type State = i64;

/// Color bitset. Bit `i` means "color i is active".
/// C++ equivalent: `unsigned long long Tag`.
pub type Tag = u64;

/// Child state vector in a transition.
pub type StateVector = Vec<State>;

// ── SymbolTag ────────────────────────────────────────────────────────────

/// A (symbol, tag) pair labeling a transition.
///
/// C++ equivalent: `SymbolTag : std::pair<Symbol, Tag>`.
/// Ordering: internal < leaf, then by symbol content, then by tag.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SymbolTag<S: SymbolTrait> {
    pub symbol: S,
    pub tag: Tag,
}

impl<S: SymbolTrait> SymbolTag<S> {
    pub fn new(symbol: S, tag: Tag) -> Self { Self { symbol, tag } }

    pub fn is_internal(&self) -> bool { self.symbol.is_internal() }
    pub fn is_leaf(&self) -> bool { self.symbol.is_leaf() }

    /// Check if a specific color bit is set.
    /// C++ equivalent: `tag(int index)`.
    pub fn has_color(&self, index: u32) -> bool {
        (self.tag >> index) & 1 != 0
    }

    /// Bitwise AND with another tag.
    /// C++ equivalent: `tag_intersection(Tag other)`.
    pub fn tag_intersection(&self, other: Tag) -> Tag {
        self.tag & other
    }
}

impl<S: SymbolTrait> Ord for SymbolTag<S> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Internal < Leaf (matching C++ operator<)
        match (self.is_internal(), other.is_internal()) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            _ => self.symbol.cmp(&other.symbol)
                     .then_with(|| self.tag.cmp(&other.tag)),
        }
    }
}

impl<S: SymbolTrait> PartialOrd for SymbolTag<S> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<S: SymbolTrait> fmt::Display for SymbolTag<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{},{}]", self.symbol, self.tag)
    }
}

// ── Transition types ─────────────────────────────────────────────────────

/// Top-down transitions: for each (symbol, tag), map parent state → set of child vectors.
///
/// C++ equivalent: `std::map<SymbolTag, std::map<State, std::set<StateVector>>>`.
pub type TopDownTransitions<S> = BTreeMap<SymbolTag<S>, BTreeMap<State, BTreeSet<StateVector>>>;

// ── Automata<S> ──────────────────────────────────────────────────────────

/// Generic non-deterministic finite tree automaton.
///
/// C++ equivalent: `AUTOQ::Automata<TT>`.
#[derive(Clone, Debug)]
pub struct Automata<S: SymbolTrait> {
    pub name: String,
    pub final_states: Vec<State>,
    pub state_num: State,
    pub qubit_num: u32,
    /// Number of symbolic variables (for SymbolicAutomata).
    pub symbolic_vars_num: u32,
    pub transitions: TopDownTransitions<S>,
    pub vars: BTreeSet<String>,
    pub constraints: String,
    pub has_loop: bool,
    pub is_topdown_deterministic: bool,
}

impl<S: SymbolTrait> Automata<S> {
    /// Create a new empty automaton for the given number of qubits.
    pub fn new(qubit_num: u32) -> Self {
        Self {
            name: String::new(),
            final_states: Vec::new(),
            state_num: 0,
            qubit_num,
            symbolic_vars_num: 0,
            transitions: BTreeMap::new(),
            vars: BTreeSet::new(),
            constraints: String::new(),
            has_loop: false,
            is_topdown_deterministic: false,
        }
    }

    // ── State management ─────────────────────────────────────────────────

    /// Allocate and return a fresh state.
    pub fn new_state(&mut self) -> State {
        let s = self.state_num;
        self.state_num += 1;
        s
    }

    // ── Transition management ────────────────────────────────────────────

    /// Add a transition: `symbol{tag}(children...) → parent`.
    pub fn add_transition(&mut self, symbol: S, tag: Tag, parent: State, children: StateVector) {
        let st = SymbolTag::new(symbol, tag);
        self.transitions
            .entry(st)
            .or_default()
            .entry(parent)
            .or_default()
            .insert(children);
    }

    // ── Query (port of query.cc) ─────────────────────────────────────────

    /// Count total number of transitions.
    ///
    /// C++ equivalent: `count_transitions()`.
    pub fn count_transitions(&self) -> usize {
        self.transitions.values()
            .flat_map(|m| m.values())
            .map(|s| s.len())
            .sum()
    }

    /// Count states actually referenced in transitions and final states.
    ///
    /// C++ equivalent: `count_states()`.
    pub fn count_states(&self) -> usize {
        let mut states = BTreeSet::new();
        for s in &self.final_states { states.insert(*s); }
        for (_, map) in &self.transitions {
            for (&parent, children_set) in map {
                states.insert(parent);
                for children in children_set {
                    for &child in children { states.insert(child); }
                }
            }
        }
        states.len()
    }

    /// Print the automaton for debugging.
    ///
    /// C++ equivalent: `print_aut(prompt)`.
    pub fn print_aut(&self, prompt: &str) {
        println!("{}Automaton '{}' (qubits={}, states={}, transitions={})",
            prompt, self.name, self.qubit_num, self.state_num, self.count_transitions());
        println!("{}Final states: {:?}", prompt, self.final_states);
        for (st, map) in &self.transitions {
            for (parent, children_set) in map {
                for children in children_set {
                    if children.is_empty() {
                        println!("{}{} → q{}", prompt, st, parent);
                    } else {
                        let ch: Vec<String> = children.iter().map(|c| format!("q{}", c)).collect();
                        println!("{}{}({}) → q{}", prompt, st, ch.join(", "), parent);
                    }
                }
            }
        }
    }
}

// ── Instance factories (port of instance.cc) ─────────────────────────────
// These create well-known automata for testing. Only `ConcreteSymbol` needs them,
// but they're generic over any S that can construct from FiveTuple-like values.

use crate::symbol::ConcreteSymbol;
use crate::fivetuple::FiveTuple;

impl Automata<ConcreteSymbol> {
    /// |0...0⟩ state: amplitude 1 on |00...0⟩, 0 elsewhere.
    ///
    /// C++ equivalent: `Automata::zero(n)` (confusing name — it means basis state |0⟩).
    pub fn zero_state(n: u32) -> Self {
        let mut aut = Self::new(n);
        let tag: Tag = 1;

        // Create leaf states
        let leaf_one = aut.new_state();  // amplitude 1 (|0⟩ path)
        let leaf_zero = aut.new_state(); // amplitude 0 (|1⟩ path)
        aut.add_transition(ConcreteSymbol::leaf(FiveTuple::one()), tag, leaf_one, vec![]);
        aut.add_transition(ConcreteSymbol::leaf(FiveTuple::zero()), tag, leaf_zero, vec![]);

        // Build tree bottom-up. For the product construction in general_single_qubit_gate
        // to work correctly, children at each qubit level must have transitions at the
        // next level (not be bare leaves). So we create a "dead" internal state at each
        // qubit level whose both children go to leaf_zero (or to the dead state at the
        // next level).
        let mut child_live = leaf_one;  // the |00...0⟩ path
        let mut child_dead = leaf_zero; // the "amplitude 0" path

        for q in (1..=n).rev() {
            let live_parent = aut.new_state();
            aut.add_transition(
                ConcreteSymbol::internal(q as i64),
                tag, live_parent,
                vec![child_live, child_dead],
            );

            // Create dead internal state: both children go to dead
            let dead_parent = aut.new_state();
            aut.add_transition(
                ConcreteSymbol::internal(q as i64),
                tag, dead_parent,
                vec![child_dead, child_dead],
            );

            child_live = live_parent;
            child_dead = dead_parent;
        }

        aut.final_states.push(child_live);
        aut
    }

    /// Uniform superposition |+⟩^⊗n: every basis state has amplitude 1/√2^n.
    ///
    /// C++ equivalent: `Automata::uniform(n)`.
    pub fn uniform(n: u32) -> Self {
        let mut aut = Self::new(n);
        let tag: Tag = 1;

        // Amplitude = 1/√2^n
        let mut amp = FiveTuple::one();
        amp.divide_by_sqrt2(n as i64);

        let leaf = aut.new_state();
        aut.add_transition(ConcreteSymbol::leaf(amp), tag, leaf, vec![]);

        // All internal nodes share the same child pair (leaf, leaf)
        let mut current = leaf;
        for q in (1..=n).rev() {
            let parent = aut.new_state();
            aut.add_transition(
                ConcreteSymbol::internal(q as i64),
                tag, parent,
                vec![current, current],
            );
            current = parent;
        }

        aut.final_states.push(current);
        aut
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_state_structure() {
        let aut = Automata::<ConcreteSymbol>::zero_state(2);
        assert_eq!(aut.qubit_num, 2);
        assert_eq!(aut.final_states.len(), 1);
        assert!(aut.count_transitions() > 0);
    }

    #[test]
    fn uniform_structure() {
        let aut = Automata::<ConcreteSymbol>::uniform(3);
        assert_eq!(aut.qubit_num, 3);
        assert_eq!(aut.final_states.len(), 1);
    }

    #[test]
    fn symbol_tag_ordering() {
        let a = SymbolTag::new(ConcreteSymbol::internal(1), 1);
        let b = SymbolTag::new(ConcreteSymbol::leaf(FiveTuple::zero()), 1);
        assert!(a < b, "internal should be less than leaf");
    }

    #[test]
    fn add_and_count_transitions() {
        let mut aut = Automata::<ConcreteSymbol>::new(1);
        let s0 = aut.new_state();
        let s1 = aut.new_state();
        aut.add_transition(ConcreteSymbol::leaf(FiveTuple::one()), 1, s1, vec![]);
        aut.add_transition(ConcreteSymbol::internal(1), 1, s0, vec![s1, s1]);
        assert_eq!(aut.count_transitions(), 2);
    }

    #[test]
    fn print_aut_smoke() {
        let aut = Automata::<ConcreteSymbol>::zero_state(1);
        aut.print_aut("test: ");
    }
}
