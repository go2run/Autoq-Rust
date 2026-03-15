# AutoQ-Rust

Rust port of [AutoQ 2.0](https://github.com/yo-han/AutoQ) — an automata-based
tool for quantum program verification using non-deterministic finite tree
automata (NFTA).

## Project Structure

```
crates/
├── autoq-core/      # Core data structures & algorithms
│   └── src/
│       ├── fivetuple.rs   # Exact complex arithmetic in ℤ[1/√2, i]
│       ├── symbol.rs      # SymbolTrait + ConcreteSymbol
│       ├── automata.rs    # Automata<S: SymbolTrait> (generic NFTA)
│       ├── gates.rs       # Quantum gate operations (X, Y, Z, H, S, T, CX, ...)
│       ├── inclusion.rs   # Language inclusion checking
│       ├── reduce.rs      # Automata reduction algorithms
│       └── general.rs     # Union, intersection, tensor product
├── autoq-parser/    # Parsing & serialization
│   └── src/
│       ├── hsl.rs         # Extended Dirac notation parser (pest)
│       └── timbuk.rs      # Timbuk format serialization
├── autoq-z3/        # Z3 SMT solver integration (optional)
│   └── src/
│       └── lib.rs         # Feature-gated Z3 bindings
└── autoq-cli/       # Command-line interface
    └── src/
        └── main.rs

tests/fixtures/      # Test data from AutoQ
├── hsl/             # HSL specification files
├── qasm/            # OpenQASM circuit files
└── lsta/            # Timbuk-format automata
```

## Building

### Basic build (no Z3)

```bash
cargo build
cargo test
```

### With Z3 support

Z3 is required for symbolic and predicate automata verification.

1. Install Z3 or locate your existing `z3.h` header and `libz3.a` library.

2. Set environment variables:
   ```bash
   export Z3_SYS_Z3_HEADER=/path/to/z3.h
   export LIBRARY_PATH=/path/to/dir/containing/libz3.a
   ```

3. Build with the feature flag:
   ```bash
   cargo build --features autoq-z3/z3-solver
   ```

### Example: using the AutoQ C++ project's Z3

If you have the original AutoQ checked out:
```bash
export Z3_SYS_Z3_HEADER=/path/to/AutoQ/include/z3/z3.h
export LIBRARY_PATH=/path/to/AutoQ
cargo build --features autoq-z3/z3-solver
```

## Architecture

### Option C+B Hybrid Type System

Inner algorithms use Rust generics (`Automata<S: SymbolTrait>`) for
zero-overhead monomorphization. The CLI layer uses enum dispatch to select
which symbol type to instantiate at runtime.

```rust
// Inner: generic, zero-cost
impl<S: SymbolTrait> Automata<S> {
    pub fn h(&mut self, t: u32) { /* Hadamard gate */ }
    pub fn cx(&mut self, c: u32, t: u32) { /* CNOT gate */ }
}

// Outer: enum dispatch at CLI
match symbol_type {
    "concrete"  => run::<ConcreteSymbol>(args),
    "symbolic"  => run::<SymbolicSymbol>(args),  // future
    "predicate" => run::<PredicateSymbol>(args),  // future
}
```

### Symbol Types

| Type | Complex Repr | Use Case | Status |
|------|-------------|----------|--------|
| `ConcreteSymbol` | `FiveTuple` | Standard verification | Implemented |
| `SymbolicSymbol` | `SymbolicComplex` | Parametric circuits | Planned |
| `PredicateSymbol` | `String` (Z3 expr) | Predicate logic | Planned |
| `IndexSymbol` | `i32` | Internal indexing | Planned |
| `ConstrainedSymbol` | `ConstrainedComplex` | Constrained verification | Planned |

## Development Notes

See [CLAUDE.md](CLAUDE.md) for detailed development notes, error logs,
and architectural decisions.
