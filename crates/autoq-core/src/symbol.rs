//! # Symbol trait and concrete implementations
//!
//! ## C++ correspondence
//! In C++ AutoQ, `Automata<TT>` is templated on `TT` which is one of:
//! - `Symbol::Concrete` (FiveTuple amplitudes)
//! - `Symbol::Symbolic` (polynomial combinations of FiveTuples)
//! - `Symbol::Predicate` (Z3 predicate strings)
//! - `Symbol::Index` (simple integer indices)
//! - `Symbol::Constrained` (constrained complex with global/local constraints)
//!
//! ## Rust design: Option C+B hybrid
//! - **Option C (generics)**: `SymbolTrait` is the trait bound on `Automata<S>`.
//!   All algorithms are monomorphized per symbol type → zero vtable overhead.
//! - **Option B (enum dispatch)**: At the CLI level, an enum selects which
//!   `Automata<ConcreteSymbol>` / `Automata<SymbolicSymbol>` / ... to use.
//!
//! Each symbol in a tree automaton is either:
//! - **Internal**: labels an internal node with a qubit index (1-based)
//! - **Leaf**: labels a leaf node with a complex amplitude

use std::fmt;
use std::hash::Hash;

/// The trait that all symbol types must implement.
///
/// This replaces C++'s template parameter `TT` in `Automata<TT>`.
/// Each symbol type defines its own complex number representation
/// (e.g., `FiveTuple` for Concrete, `SymbolicComplex` for Symbolic).
///
/// ## Required operations
/// Gate algorithms need arithmetic on leaf amplitudes (`add`, `mul`, `negate`,
/// `divide_by_sqrt2`, `counterclockwise`). These are defined on the associated
/// `Complex` type.
pub trait SymbolTrait: Clone + Eq + Ord + Hash + fmt::Debug + fmt::Display + Send + Sync {
    /// The amplitude representation used by leaf nodes.
    /// - `FiveTuple` for Concrete
    /// - `SymbolicComplex` for Symbolic (future)
    /// - `String` for Predicate (future)
    type Complex: Clone + Eq + Ord + Hash + fmt::Debug + fmt::Display + Send + Sync;

    // --- Construction ---

    /// Create an internal node symbol for the given qubit index (1-based).
    ///
    /// C++ equivalent: `Concrete(qubit)` with `internal = true`.
    fn internal(qubit: i64) -> Self;

    /// Create a leaf node symbol with the given amplitude.
    ///
    /// C++ equivalent: `Concrete(complex)` with `internal = false`.
    fn leaf(complex: Self::Complex) -> Self;

    // --- Classification ---

    fn is_internal(&self) -> bool;
    fn is_leaf(&self) -> bool { !self.is_internal() }

    /// Get the qubit index. Panics if called on a leaf.
    ///
    /// C++ equivalent: `symbol.qubit()`.
    fn qubit(&self) -> i64;

    /// Get a reference to the leaf amplitude. Panics if called on an internal node.
    fn complex(&self) -> &Self::Complex;

    /// Consume and return the leaf amplitude. Panics if called on an internal node.
    fn into_complex(self) -> Self::Complex;

    // --- Amplitude arithmetic (delegated to Complex) ---
    // These exist on the trait so gate algorithms can operate generically.

    fn complex_zero() -> Self::Complex;
    fn complex_one() -> Self::Complex;
    fn complex_add(a: &Self::Complex, b: &Self::Complex) -> Self::Complex;
    fn complex_sub(a: &Self::Complex, b: &Self::Complex) -> Self::Complex;
    fn complex_mul(a: &Self::Complex, b: &Self::Complex) -> Self::Complex;
    fn complex_neg(c: &Self::Complex) -> Self::Complex;
    fn complex_is_zero(c: &Self::Complex) -> bool;

    /// Divide by √2^times. Used by Hadamard and other gates.
    ///
    /// C++ equivalent: `complex.divide_by_the_square_root_of_two(times)`.
    fn complex_divide_by_sqrt2(c: &mut Self::Complex, times: i64);

    /// Rotate counterclockwise by θ = num/den turns (multiples of full circle).
    /// Requires 8·θ ∈ ℤ (Clifford+T angles only).
    ///
    /// C++ equivalent: `complex.counterclockwise(theta)`.
    fn complex_counterclockwise(c: &mut Self::Complex, theta_num: i64, theta_den: i64);

    /// Rotate clockwise by θ = num/den turns.
    fn complex_clockwise(c: &mut Self::Complex, theta_num: i64, theta_den: i64);

    /// Multiply by cos(π·θ).
    /// C++ equivalent: `complex.multiply_cos(theta)`.
    fn complex_multiply_cos(c: &Self::Complex, theta_num: i64, theta_den: i64) -> Self::Complex;

    /// Multiply by i·sin(π·θ).
    fn complex_multiply_isin(c: &Self::Complex, theta_num: i64, theta_den: i64) -> Self::Complex;

    /// Simplify the fraction representation if applicable.
    /// No-op for symbol types that don't support it (e.g., Predicate).
    fn complex_fraction_simplification(c: &mut Self::Complex) { let _ = c; }

    /// Set amplitude to zero. Used by `back_to_zero()` in C++.
    fn complex_back_to_zero(c: &mut Self::Complex);
}

// ═══════════════════════════════════════════════════════════════════════════
// ConcreteSymbol: the primary symbol type using FiveTuple
// ═══════════════════════════════════════════════════════════════════════════

use crate::fivetuple::FiveTuple;

/// Concrete symbol: internal nodes carry qubit indices, leaves carry FiveTuple amplitudes.
///
/// C++ equivalent: `AUTOQ::Symbol::Concrete` in `include/autoq/symbol/concrete.hh`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ConcreteSymbol {
    Internal(i64),
    Leaf(FiveTuple),
}

/// Ordering: internal < leaf (matching C++ `operator<`), then by content.
impl Ord for ConcreteSymbol {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (ConcreteSymbol::Internal(_), ConcreteSymbol::Leaf(_)) => std::cmp::Ordering::Less,
            (ConcreteSymbol::Leaf(_), ConcreteSymbol::Internal(_)) => std::cmp::Ordering::Greater,
            (ConcreteSymbol::Internal(a), ConcreteSymbol::Internal(b)) => a.cmp(b),
            (ConcreteSymbol::Leaf(a), ConcreteSymbol::Leaf(b)) => a.cmp(b),
        }
    }
}

impl PartialOrd for ConcreteSymbol {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for ConcreteSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConcreteSymbol::Internal(q) => write!(f, "x{}", q),
            ConcreteSymbol::Leaf(ft) => write!(f, "{}", ft),
        }
    }
}

impl SymbolTrait for ConcreteSymbol {
    type Complex = FiveTuple;

    fn internal(qubit: i64) -> Self { ConcreteSymbol::Internal(qubit) }
    fn leaf(complex: FiveTuple) -> Self { ConcreteSymbol::Leaf(complex) }

    fn is_internal(&self) -> bool {
        matches!(self, ConcreteSymbol::Internal(_))
    }

    fn qubit(&self) -> i64 {
        match self {
            ConcreteSymbol::Internal(q) => *q,
            ConcreteSymbol::Leaf(_) => panic!("Leaf symbols do not have qubit()"),
        }
    }

    fn complex(&self) -> &FiveTuple {
        match self {
            ConcreteSymbol::Leaf(ft) => ft,
            ConcreteSymbol::Internal(_) => panic!("Internal symbols do not have complex()"),
        }
    }

    fn into_complex(self) -> FiveTuple {
        match self {
            ConcreteSymbol::Leaf(ft) => ft,
            ConcreteSymbol::Internal(_) => panic!("Internal symbols do not have complex()"),
        }
    }

    fn complex_zero() -> FiveTuple { FiveTuple::zero() }
    fn complex_one() -> FiveTuple { FiveTuple::one() }

    fn complex_add(a: &FiveTuple, b: &FiveTuple) -> FiveTuple { a.clone() + b.clone() }
    fn complex_sub(a: &FiveTuple, b: &FiveTuple) -> FiveTuple { a.clone() - b.clone() }
    fn complex_mul(a: &FiveTuple, b: &FiveTuple) -> FiveTuple { a.clone() * b.clone() }
    fn complex_neg(c: &FiveTuple) -> FiveTuple { -c.clone() }
    fn complex_is_zero(c: &FiveTuple) -> bool { c.is_zero() }

    fn complex_divide_by_sqrt2(c: &mut FiveTuple, times: i64) {
        c.divide_by_sqrt2(times);
    }

    fn complex_counterclockwise(c: &mut FiveTuple, theta_num: i64, theta_den: i64) {
        c.counterclockwise(theta_num, theta_den);
    }

    fn complex_clockwise(c: &mut FiveTuple, theta_num: i64, theta_den: i64) {
        c.clockwise(theta_num, theta_den);
    }

    fn complex_multiply_cos(c: &FiveTuple, theta_num: i64, theta_den: i64) -> FiveTuple {
        c.multiply_cos(theta_num, theta_den)
    }

    fn complex_multiply_isin(c: &FiveTuple, theta_num: i64, theta_den: i64) -> FiveTuple {
        c.multiply_isin(theta_num, theta_den)
    }

    fn complex_fraction_simplification(c: &mut FiveTuple) {
        c.fraction_simplification();
    }

    fn complex_back_to_zero(c: &mut FiveTuple) {
        c.back_to_zero();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concrete_internal_ordering() {
        let q1 = ConcreteSymbol::internal(1);
        let q2 = ConcreteSymbol::internal(2);
        let leaf = ConcreteSymbol::leaf(FiveTuple::zero());
        // Internal < Leaf
        assert!(q1 < leaf);
        // Internal ordering by qubit index
        assert!(q1 < q2);
    }

    #[test]
    fn concrete_leaf_operations() {
        let a = FiveTuple::one();
        let b = FiveTuple::inv_sqrt2();
        let sum = ConcreteSymbol::complex_add(&a, &b);
        assert!(!sum.is_zero());
    }

    #[test]
    fn concrete_symbol_display() {
        assert_eq!(format!("{}", ConcreteSymbol::internal(3)), "x3");
    }

    #[test]
    fn concrete_roundtrip() {
        let ft = FiveTuple::inv_sqrt2();
        let sym = ConcreteSymbol::leaf(ft.clone());
        assert!(sym.is_leaf());
        assert_eq!(sym.complex(), &ft);
        assert_eq!(sym.into_complex(), ft);
    }
}
