# CLAUDE.md — Development Notes for AutoQ-Rust

## Project Overview
AutoQ-Rust is a faithful Rust port of AutoQ 2.0, a C++ automata-based tool for
quantum program verification using non-deterministic finite tree automata (NFTA).

Original C++ source: `/home/user/AutoQ/`
Prior Rust prototype: `/home/user/AutoQ/rust_proto/` (8-crate workspace)

## Architecture Decisions

### Option C+B Hybrid (Symbol type system)
- **Decision**: Use generics for inner algorithms + enum dispatch at CLI
- **Why**: Zero vtable overhead in hot paths; compiler optimizes monomorphized code
- **Trade-off**: Binary size grows per symbol type; compile time increases
- **Alternative rejected**: Pure enum dispatch (Option B only) — would add match
  overhead in every gate/inclusion operation

### 4-Crate Structure (not 8)
- `autoq-core`: all core types + algorithms in one crate
- `autoq-parser`: all parsers (HSL, QASM, Timbuk)
- `autoq-z3`: Z3 integration, feature-gated
- `autoq-cli`: binary entry point
- **Why**: The original 8-crate structure forced orphan-rule workarounds
  (extension traits for gate operations). With 4 crates, gates can `impl Automata<S>`
  directly.

### Z3 Feature Gate
- `autoq-z3` compiles without Z3 by default
- Enable with `--features z3-solver`
- **Why**: Z3 is a 50MB dependency with complex build requirements.
  Concrete symbol operations (the common case) don't need it.

## Errors & Attempts Log

### FiveTuple Division
- **Attempt 1**: Tried using `num-rational::BigRational` for division.
  Failed because FiveTuple is not a rational number — it's an element of ℤ[1/√2, i].
- **Solution**: Ported the C++ conjugate-multiplication algorithm directly.
  Key insight: multiply by the algebraic conjugate to make denominator a rational
  integer, then divide out GCD and extract powers of 2.
- **Gotcha**: BigInt type inference fails with mixed `2*&a*&b` expressions.
  Must annotate intermediate types explicitly: `BigInt::from(2)*&a*&b`.

### Symbol Trait Design
- **Attempt 1**: Tried putting arithmetic methods directly on `SymbolTrait`.
  This works but makes the trait very large (15+ methods).
- **Decision**: Keep arithmetic on the trait as associated methods (static dispatch).
  Alternative would be a separate `ComplexArithmetic` trait, but that adds complexity
  for minimal benefit.
- **Note**: `complex_multiply_cos` and `complex_multiply_isin` take `(i64, i64)`
  instead of `boost::rational<cpp_int>`. This is sufficient for Clifford+T gates
  (all angles are multiples of π/4). For arbitrary rotation gates, we may need
  `num_rational::BigRational` in the future.

### Automata Generics
- **Attempt 1**: Tried `Automata<S>` where `S` is just the leaf type.
  Problem: internal nodes also need to be `S` for uniform transition maps.
- **Solution**: `SymbolTrait` has both `internal(qubit)` and `leaf(complex)` constructors.
  `ConcreteSymbol` is an enum with `Internal(i64)` and `Leaf(FiveTuple)` variants.

### CZ Gate
- **Attempt 1**: Tried to use `states_on_one_branch` intersection for CZ.
  Problem: intersection doesn't correctly identify states reachable through
  BOTH control=|1⟩ AND target=|1⟩ because the automaton is a tree, not a DAG.
- **Current**: Simplified approach checks if leaf state appears on both branches.
  This works for path-structured automata (deterministic per color).
  May need revision for highly nondeterministic cases.

### Pest vs ANTLR4
- ANTLR4's semantic predicates (`{isNonZero($N.text)}?`) have no pest equivalent.
- Solution: parse permissively, validate in Rust code after parsing.
- The symbolic ket format (`|ss0001>`) and sum quantifiers (`∑`) are not yet
  covered by the pest grammar. These require significant grammar extensions.

## C++ to Rust Mapping

| C++ | Rust | Notes |
|-----|------|-------|
| `Automata<TT>` | `Automata<S: SymbolTrait>` | Generic over symbol type |
| `Symbol::Concrete` | `ConcreteSymbol` | Enum: Internal(i64) / Leaf(FiveTuple) |
| `Complex::FiveTuple` | `FiveTuple` | Named fields instead of inheriting Vec |
| `boost::multiprecision::cpp_int` | `num_bigint::BigInt` | |
| `boost::rational<cpp_int>` | `num_rational::BigRational` | Used in rotation angles |
| `std::map<K,V>` | `BTreeMap<K,V>` | Deterministic ordering |
| `unsigned long long Tag` | `u64 Tag` | Color bitset |
| `ANTLR4 runtime` | `pest` PEG parser | |
| `CLI11` | `clap` | |
| `inline static` globals | Explicit parameters or `thread_local!` | No hidden global state |

## Test Data
Test fixtures in `tests/fixtures/` are copied from AutoQ's unit test cases:
- `hsl/`: HSL specification files (pre/post conditions)
- `qasm/`: OpenQASM circuit files
- `lsta/`: Timbuk-format automata serializations

## Next Steps (Priority Order)
1. Port HSL parser from prototype (adapt for generic Automata<S>)
2. Port Timbuk serialization/deserialization
3. Implement QASM parser with pest
4. Implement Rx(θ), Ry(θ), Rz(θ), CCX gates
5. Implement `sim_reduce` (requires LTS simulation algorithm)
6. Implement `execute()` (QASM circuit execution engine)
7. Z3 integration for Symbolic/Predicate symbol types
8. Loop summarization
