//! # autoq-z3 — Z3 SMT solver integration
//!
//! ## C++ correspondence
//! Z3 is used in AutoQ for:
//! 1. Symbolic amplitude checking (SymbolicAutomata inclusion)
//! 2. Predicate automata verification
//! 3. Loop summarization convergence conditions
//! 4. FiveTuple real/imag SMT expression generation
//!
//! ## Feature gating
//! This crate is **optional**. The `z3-solver` feature enables Z3 bindings.
//! Without it, only stub types are available.
//!
//! ```toml
//! # Enable Z3:
//! autoq-z3 = { path = "crates/autoq-z3", features = ["z3-solver"] }
//! ```
//!
//! ## Build requirements (when z3-solver enabled)
//! Set environment variables before building:
//! ```bash
//! export Z3_SYS_Z3_HEADER=/path/to/z3.h
//! export LIBRARY_PATH=/path/to/dir/containing/libz3.a
//! cargo build --features z3-solver
//! ```

#[cfg(feature = "z3-solver")]
pub mod solver;

#[cfg(not(feature = "z3-solver"))]
pub mod solver {
    //! Stub module when Z3 is not available.
    //! Provides type signatures but panics at runtime.

    /// Placeholder for when Z3 is not compiled in.
    pub struct SatChecker;

    impl SatChecker {
        pub fn new() -> Self {
            panic!("Z3 support not compiled. Rebuild with: cargo build --features z3-solver")
        }
    }
}
