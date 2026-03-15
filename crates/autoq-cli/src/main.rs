//! # autoq — CLI for quantum program verification
//!
//! ## C++ correspondence
//! Port of `cli/autoq.cc` (400+ lines).
//!
//! ## Commands
//! - `ex`: execute a circuit with precondition
//! - `ver`: verify circuit against pre/post conditions
//! - `eq`: check equivalence of two circuits
//! - `print`: print a quantum state set from HSL
//!
//! ## Status
//! Stub — skeleton only.

fn main() {
    println!("autoq-rust: not yet implemented");
    println!("Available crates:");
    println!("  autoq-core   — core automata, gates, inclusion, reduction");
    println!("  autoq-parser — HSL, Timbuk, QASM parsers");
    println!("  autoq-z3     — Z3 SMT solver integration (feature-gated)");
}
