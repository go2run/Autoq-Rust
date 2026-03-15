//! # autoq-core — Core data structures and algorithms for AutoQ
//!
//! This crate contains:
//! - `FiveTuple`: exact complex arithmetic in ℤ[1/√2, i]
//! - `SymbolTrait`: the trait that all symbol types must implement
//! - `Automata<S>`: generic non-deterministic finite tree automaton
//! - Gate operations, inclusion checking, reduction algorithms
//!
//! ## Architecture: Option C+B hybrid
//! Inner algorithms use generics (`Automata<S: SymbolTrait>`) for zero-overhead
//! monomorphization. The CLI layer uses enum dispatch to select the symbol type
//! at runtime.
//!
//! ## Correspondence to C++ AutoQ
//! - `fivetuple` → `include/autoq/complex/fivetuple.hh`
//! - `symbol` → `include/autoq/symbol/*.hh`
//! - `automata` → `include/autoq/aut_description.hh`
//! - `gates` → `src/gate.cc`
//! - `inclusion` → `src/inclusion.cc`
//! - `reduce` → `src/reduce.cc`
//! - `general` → `src/general.cc`

pub mod fivetuple;
pub mod symbol;
pub mod automata;
pub mod gates;
pub mod inclusion;
pub mod reduce;
pub mod general;

pub use fivetuple::FiveTuple;
pub use symbol::SymbolTrait;
pub use automata::Automata;
