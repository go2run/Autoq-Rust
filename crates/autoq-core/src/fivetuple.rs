//! # FiveTuple — Exact complex arithmetic in ℤ[1/√2, i]
//!
//! ## C++ correspondence
//! Faithful port of `include/autoq/complex/fivetuple.hh`.
//!
//! ## Representation
//! `FiveTuple { a, b, c, d, k }` encodes the complex number:
//! ```text
//! (a + b·(1+i)/√2 + c·i + d·(-1+i)/√2) / √2^k
//! ```
//! The four basis elements {1, (1+i)/√2, i, (-1+i)/√2} correspond to the
//! eight 8th-roots of unity. Every Clifford+T amplitude lives in this ring.
//!
//! ## Differences from C++
//! - C++ inherits `std::vector<cpp_int>`; Rust uses named fields for clarity.
//! - `Ord` is implemented for use as `BTreeMap` keys (C++ uses `operator<`).
//! - `Div` is implemented here (C++ also has `operator/` in fivetuple.hh).
//! - Rotation methods take `(i64, i64)` instead of `boost::rational`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Exact complex number in ℤ[1/√2, i].
///
/// Invariant: value = (a + b·ω + c·ω² + d·ω³) / √2^k
/// where ω = (1+i)/√2 = e^{iπ/4}.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct FiveTuple {
    pub a: BigInt,
    pub b: BigInt,
    pub c: BigInt,
    pub d: BigInt,
    /// Exponent of √2 in the denominator. May be negative.
    pub k: BigInt,
}

// ── Ordering (for BTreeMap keys) ─────────────────────────────────────────
// C++ orders by vector comparison. We replicate: compare k first, then a,b,c,d.
// This is NOT mathematical ordering — it's a canonical key ordering.

impl Ord for FiveTuple {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.k.cmp(&other.k)
            .then_with(|| self.a.cmp(&other.a))
            .then_with(|| self.b.cmp(&other.b))
            .then_with(|| self.c.cmp(&other.c))
            .then_with(|| self.d.cmp(&other.d))
    }
}

impl PartialOrd for FiveTuple {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// ── Constants ────────────────────────────────────────────────────────────

impl FiveTuple {
    pub fn zero() -> Self {
        Self { a: BigInt::zero(), b: BigInt::zero(), c: BigInt::zero(), d: BigInt::zero(), k: BigInt::zero() }
    }

    pub fn one() -> Self {
        Self { a: BigInt::one(), b: BigInt::zero(), c: BigInt::zero(), d: BigInt::zero(), k: BigInt::zero() }
    }

    /// The imaginary unit i = [0,0,1,0,0].
    pub fn imag_unit() -> Self {
        Self { a: BigInt::zero(), b: BigInt::zero(), c: BigInt::one(), d: BigInt::zero(), k: BigInt::zero() }
    }

    /// 1/√2 = [1,0,0,0,1].
    pub fn inv_sqrt2() -> Self {
        Self { a: BigInt::one(), b: BigInt::zero(), c: BigInt::zero(), d: BigInt::zero(), k: BigInt::one() }
    }

    /// √2 = [1,0,0,0,-1].
    ///
    /// C++ equivalent: `FiveTuple(1).divide_by_the_square_root_of_two(-1)`.
    pub fn sqrt2() -> Self {
        Self { a: BigInt::one(), b: BigInt::zero(), c: BigInt::zero(), d: BigInt::zero(), k: -BigInt::one() }
    }

    pub fn from_int(n: impl Into<BigInt>) -> Self {
        Self { a: n.into(), b: BigInt::zero(), c: BigInt::zero(), d: BigInt::zero(), k: BigInt::zero() }
    }

    /// e^{iπ·θ} where θ = num/den.
    pub fn angle(theta_num: i64, theta_den: i64) -> Self {
        let mut v = Self::one();
        v.counterclockwise(theta_num, theta_den);
        v
    }

    pub fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero() && self.c.is_zero() && self.d.is_zero()
    }

    /// Set all coefficients to zero, keeping k unchanged.
    /// C++ equivalent: `back_to_zero()`.
    pub fn back_to_zero(&mut self) {
        self.a = BigInt::zero();
        self.b = BigInt::zero();
        self.c = BigInt::zero();
        self.d = BigInt::zero();
    }
}

// ── k-manipulation ───────────────────────────────────────────────────────

impl FiveTuple {
    /// Multiply numerator by √2 and increment k by 1, preserving value.
    ///
    /// C++ equivalent: `increase_k()`.
    /// Transform: [a,b,c,d,k] → [b-d, a+c, b+d, c-a, k+1]
    pub fn increase_k(&mut self) {
        let (a, b, c, d) = (self.a.clone(), self.b.clone(), self.c.clone(), self.d.clone());
        self.a = &b - &d;
        self.b = &a + &c;
        self.c = &b + &d;
        self.d = &c - &a;
        self.k += 1;
    }

    /// Increase k until it reaches target.
    pub fn increase_to_k(&mut self, target: &BigInt) {
        while &self.k < target {
            self.increase_k();
        }
    }

    /// Divide value by √2^times (increment k).
    ///
    /// C++ equivalent: `divide_by_the_square_root_of_two(times)`.
    pub fn divide_by_sqrt2(&mut self, times: i64) {
        self.k += times;
    }

    /// Cancel common factors of 2: while a,b,c,d all even and k≥2, halve and subtract 2 from k.
    ///
    /// C++ equivalent: `fraction_simplification()`.
    pub fn fraction_simplification(&mut self) {
        if self.is_zero() {
            self.k = BigInt::zero();
            return;
        }
        let two = BigInt::from(2);
        while self.k >= two {
            if self.a.is_even() && self.b.is_even() && self.c.is_even() && self.d.is_even() {
                self.a >>= 1;
                self.b >>= 1;
                self.c >>= 1;
                self.d >>= 1;
                self.k -= 2;
            } else {
                break;
            }
        }
    }
}

// ── Rotation ─────────────────────────────────────────────────────────────

impl FiveTuple {
    /// Rotate counterclockwise by θ = num/den turns.
    /// One CCW step: [a,b,c,d] → [-d, a, b, c].
    ///
    /// C++ equivalent: `counterclockwise(theta)`.
    pub fn counterclockwise(&mut self, theta_num: i64, theta_den: i64) {
        let r8 = 8i64 * theta_num;
        assert!(r8 % theta_den == 0, "8·θ must be integer, got {}/{}", theta_num, theta_den);
        let mut r = (r8 / theta_den).rem_euclid(8);
        while r > 0 {
            let tmp = self.d.clone();
            self.d = self.c.clone();
            self.c = self.b.clone();
            self.b = self.a.clone();
            self.a = -tmp;
            r -= 1;
        }
    }

    /// Rotate clockwise by θ = num/den turns.
    /// One CW step: [a,b,c,d] → [b, c, d, -a].
    ///
    /// C++ equivalent: `clockwise(theta)`.
    pub fn clockwise(&mut self, theta_num: i64, theta_den: i64) {
        let r8 = 8i64 * theta_num;
        assert!(r8 % theta_den == 0, "8·θ must be integer");
        let mut r = (r8 / theta_den).rem_euclid(8);
        while r > 0 {
            let tmp = self.a.clone();
            self.a = self.b.clone();
            self.b = self.c.clone();
            self.c = self.d.clone();
            self.d = -tmp;
            r -= 1;
        }
    }

    /// self · cos(π·θ) via Euler: cos = (e^{iθ} + e^{-iθ}) / 2.
    ///
    /// C++ equivalent: `multiply_cos(theta)`.
    pub fn multiply_cos(&self, theta_num: i64, theta_den: i64) -> Self {
        let mut c1 = self.clone();
        let mut c2 = self.clone();
        c1.counterclockwise(theta_num, theta_den);
        c2.clockwise(theta_num, theta_den);
        let mut result = c1 + c2;
        result.divide_by_sqrt2(2);
        result
    }

    /// self · i·sin(π·θ) via Euler: i·sin = (e^{iθ} - e^{-iθ}) / 2.
    ///
    /// C++ equivalent: `multiply_isin(theta)`.
    pub fn multiply_isin(&self, theta_num: i64, theta_den: i64) -> Self {
        let mut c1 = self.clone();
        let mut c2 = self.clone();
        c1.counterclockwise(theta_num, theta_den);
        c2.clockwise(theta_num, theta_den);
        let mut result = c1 - c2;
        result.divide_by_sqrt2(2);
        result
    }
}

// ── Arithmetic ───────────────────────────────────────────────────────────

/// Internal: align two FiveTuples to the same k, then add or subtract.
fn binary_op(mut lhs: FiveTuple, mut rhs: FiveTuple, add: bool) -> FiveTuple {
    while lhs.k < rhs.k { lhs.increase_k(); }
    while rhs.k < lhs.k { rhs.increase_k(); }
    let k = lhs.k.clone();
    FiveTuple {
        a: if add { lhs.a + rhs.a } else { lhs.a - rhs.a },
        b: if add { lhs.b + rhs.b } else { lhs.b - rhs.b },
        c: if add { lhs.c + rhs.c } else { lhs.c - rhs.c },
        d: if add { lhs.d + rhs.d } else { lhs.d - rhs.d },
        k,
    }
}

impl Add for FiveTuple {
    type Output = Self;
    fn add(self, rhs: Self) -> Self { binary_op(self, rhs, true) }
}

impl Sub for FiveTuple {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self { binary_op(self, rhs, false) }
}

impl Neg for FiveTuple {
    type Output = Self;
    fn neg(self) -> Self {
        FiveTuple { a: -self.a, b: -self.b, c: -self.c, d: -self.d, k: self.k }
    }
}

/// Multiplication in ℤ[ω] where ω⁴ = -1.
///
/// C++ equivalent: `operator*`.
impl Mul for FiveTuple {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        FiveTuple {
            a: &self.a * &r.a - &self.b * &r.d - &self.c * &r.c - &self.d * &r.b,
            b: &self.a * &r.b + &self.b * &r.a - &self.c * &r.d - &self.d * &r.c,
            c: &self.a * &r.c + &self.b * &r.b + &self.c * &r.a - &self.d * &r.d,
            d: &self.a * &r.d + &self.b * &r.c + &self.c * &r.b + &self.d * &r.a,
            k: &self.k + &r.k,
        }
    }
}

/// Division in ℤ[1/√2, i].
///
/// C++ equivalent: `operator/` in fivetuple.hh.
/// Uses the conjugate-multiplication trick: multiply numerator and denominator
/// by the "conjugate" to make the denominator a rational integer, then simplify.
///
/// The conjugate of [a,b,c,d] is [a³+2abd+ac²+cd²-b²c, ..., ..., ...].
impl Div for FiveTuple {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let (a, b, c, d) = (rhs.a.clone(), rhs.b.clone(), rhs.c.clone(), rhs.d.clone());

        // Build conjugate of rhs
        let mut conj = FiveTuple::zero();
        conj.a = &a*&a*&a + 2*&a*&b*&d + &a*&c*&c + &c*&d*&d - &b*&b*&c;
        conj.b = -(&a*&a*&b) - &b*&b*&d + &b*&c*&c - &d*&d*&d - 2*&a*&c*&d;
        conj.c = -(&a*&a*&c) + 2*&b*&c*&d - &c*&c*&c - &a*&d*&d + &a*&b*&b;
        conj.d = -(&b*&d*&d) + 2*&a*&b*&c - &b*&b*&b - &a*&a*&d + &c*&c*&d;
        conj.k = -rhs.k;

        let mut result = self * conj;

        // Denominator after conjugation: a real non-negative integer times a power of 2
        let term1: BigInt = &a*&a + BigInt::from(2)*&b*&d - &c*&c;
        let term2: BigInt = &d*&d + BigInt::from(2)*&a*&c - &b*&b;
        let mut dd: BigInt = (&term1 * &term1) + (&term2 * &term2);

        // GCD simplification
        let cf: BigInt = [&result.a, &result.b, &result.c, &result.d].iter()
            .fold(dd.clone(), |acc: BigInt, x| acc.gcd(x));
        result.a /= &cf;
        result.b /= &cf;
        result.c /= &cf;
        result.d /= &cf;
        dd /= &cf;

        // Extract powers of 2 from dd into k
        while (&dd & BigInt::one()).is_zero() && !dd.is_zero() {
            result.k += 2;
            dd >>= 1;
        }
        assert!(dd.abs() == BigInt::one(), "Denominator must reduce to ±1, got {}", dd);

        if dd.is_negative() {
            result.a = -result.a;
            result.b = -result.b;
            result.c = -result.c;
            result.d = -result.d;
        }

        result
    }
}

// ── Value equality ───────────────────────────────────────────────────────

impl FiveTuple {
    /// Compare values after normalising to the same k.
    /// Unlike `==` which compares representation, this checks mathematical equality.
    ///
    /// C++ equivalent: `valueEqual(o)`.
    pub fn value_eq(&self, other: &Self) -> bool {
        let mut lhs = self.clone();
        let mut rhs = other.clone();
        while lhs.k < rhs.k { lhs.increase_k(); }
        while rhs.k < lhs.k { rhs.increase_k(); }
        lhs.a == rhs.a && lhs.b == rhs.b && lhs.c == rhs.c && lhs.d == rhs.d
    }

    /// Extract real part as a FiveTuple: {b-d, a, 0, -a, k+1}.
    ///
    /// C++ equivalent: `real()`.
    pub fn real_part(&self) -> Self {
        FiveTuple {
            a: &self.b - &self.d, b: self.a.clone(),
            c: BigInt::zero(), d: -self.a.clone(), k: &self.k + 1,
        }
    }

    /// Extract imaginary part: {0, c, b+d, c, k+1}.
    ///
    /// C++ equivalent: `imag()`.
    pub fn imag_part(&self) -> Self {
        FiveTuple {
            a: BigInt::zero(), b: self.c.clone(),
            c: &self.b + &self.d, d: self.c.clone(), k: &self.k + 1,
        }
    }

    /// Convert to rational if purely real and rational.
    ///
    /// C++ equivalent: `to_rational()`.
    pub fn to_rational(&self) -> Option<num_rational::BigRational> {
        if !self.imag_part().is_zero() { return None; }
        let real = self.real_part();
        let k_i64 = real.k.to_i64()?;

        // Real part = (b-d)/√2^(k+1) + a/√2^k
        // For this to be rational, irrational terms must cancel
        let b_minus_d = &real.a; // after real_part transform, a = orig.b - orig.d
        let a_coeff = &real.b;   // after real_part transform, b = orig.a

        if (k_i64 + 1) % 2 != 0 && !b_minus_d.is_zero() { return None; }
        if k_i64 % 2 != 0 && !a_coeff.is_zero() { return None; }

        use num_rational::BigRational;
        let two = BigInt::from(2);

        if k_i64 >= 0 {
            let term1 = if !b_minus_d.is_zero() {
                BigRational::new(b_minus_d.clone(), two.pow(((k_i64 + 1) / 2) as u32))
            } else { BigRational::zero() };

            let term2 = if !a_coeff.is_zero() {
                BigRational::new(a_coeff.clone(), two.pow((k_i64 / 2) as u32))
            } else { BigRational::zero() };

            Some(term1 + term2)
        } else {
            let abs_k = (-k_i64) as u32;
            let term1 = BigRational::from(b_minus_d.clone() * two.pow((abs_k.saturating_sub(1)) / 2));
            let term2 = BigRational::from(a_coeff.clone() * two.pow(abs_k / 2));
            Some(term1 + term2)
        }
    }
}

// ── Display ──────────────────────────────────────────────────────────────

impl fmt::Display for FiveTuple {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() { return write!(f, "0"); }
        write!(f, "[{},{},{},{},{}]", self.a, self.b, self.c, self.d, self.k)
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(n: i64) -> BigInt { BigInt::from(n) }

    #[test]
    fn zero_is_zero() { assert!(FiveTuple::zero().is_zero()); }

    #[test]
    fn one_plus_one() {
        assert!((FiveTuple::one() + FiveTuple::one()).value_eq(&FiveTuple::from_int(2)));
    }

    #[test]
    fn inv_sqrt2_squared_is_half() {
        let s = FiveTuple::inv_sqrt2();
        let half = FiveTuple { a: bi(1), b: bi(0), c: bi(0), d: bi(0), k: bi(2) };
        assert!((s.clone() * s).value_eq(&half));
    }

    #[test]
    fn multiply_by_i() {
        let v = FiveTuple { a: bi(3), b: bi(1), c: bi(2), d: bi(5), k: bi(0) };
        let result = v * FiveTuple::imag_unit();
        let expected = FiveTuple { a: bi(-2), b: bi(-5), c: bi(3), d: bi(1), k: bi(0) };
        assert!(result.value_eq(&expected));
    }

    #[test]
    fn ccw_full_turn_identity() {
        let v = FiveTuple { a: bi(3), b: bi(1), c: bi(2), d: bi(5), k: bi(0) };
        let mut r = v.clone();
        r.counterclockwise(1, 1);
        assert!(r.value_eq(&v));
    }

    #[test]
    fn increase_k_preserves_value() {
        let v = FiveTuple::one();
        let mut u = v.clone();
        u.increase_k();
        assert!(u.value_eq(&v));
    }

    #[test]
    fn fraction_simplification_works() {
        let mut v = FiveTuple { a: bi(4), b: bi(4), c: bi(4), d: bi(4), k: bi(4) };
        let orig = v.clone();
        v.fraction_simplification();
        assert!(v.value_eq(&orig));
        assert_eq!(v.k, bi(0));
        assert_eq!(v.a, bi(1));
    }

    #[test]
    fn division_basic() {
        // (1/√2) / (1/√2) = 1
        let a = FiveTuple::inv_sqrt2();
        let result = a.clone() / a;
        assert!(result.value_eq(&FiveTuple::one()), "got {:?}", result);
    }

    #[test]
    fn division_by_one() {
        let v = FiveTuple { a: bi(3), b: bi(1), c: bi(2), d: bi(5), k: bi(0) };
        let result = v.clone() / FiveTuple::one();
        assert!(result.value_eq(&v), "got {:?}", result);
    }

    #[test]
    fn division_complex() {
        // i / i = 1
        let i = FiveTuple::imag_unit();
        let result = i.clone() / i;
        assert!(result.value_eq(&FiveTuple::one()), "i/i should be 1, got {:?}", result);
    }

    #[test]
    fn multiply_cos_pi_over_4() {
        let result = FiveTuple::one().multiply_cos(1, 8);
        assert!(result.value_eq(&FiveTuple::inv_sqrt2()));
    }

    #[test]
    fn ordering_deterministic() {
        let a = FiveTuple::zero();
        let b = FiveTuple::one();
        assert!(a < b || b < a); // Just ensure Ord works
    }

    #[test]
    fn t_gate_amplitude() {
        let t = FiveTuple::angle(1, 8);
        let expected = FiveTuple { a: bi(0), b: bi(1), c: bi(0), d: bi(0), k: bi(0) };
        assert!(t.value_eq(&expected));
    }
}
