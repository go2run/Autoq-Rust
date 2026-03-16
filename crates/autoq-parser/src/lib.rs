//! # autoq-parser — Parsing and serialization for AutoQ
//!
//! ## Modules
//! - `hsl`: Extended Dirac notation parser (replaces ANTLR4)
//! - `timbuk`: Timbuk format serialization/deserialization
//! - `qasm`: OpenQASM 2.0 line-by-line gate executor

pub mod hsl;
pub mod timbuk;
pub mod qasm;
