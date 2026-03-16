//! # autoq — CLI for quantum program verification
//!
//! ## C++ correspondence
//! Port of `cli/autoq.cc`.
//!
//! ## Commands
//! - `ex`:    execute a circuit with precondition
//! - `ver`:   verify circuit against pre/post conditions
//! - `eq`:    check equivalence of two circuits
//! - `print`: print a quantum state set from HSL

use anyhow::{Context, Result};
use autoq_core::inclusion;
use autoq_core::symbol::ConcreteSymbol;
use autoq_core::Automata;
use autoq_parser::{hsl, qasm, timbuk};
use clap::{Parser, Subcommand};
use std::fs;
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "autoq",
    about = "AutoQ 2.0: An automata-based tool for quantum program verification.",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Execute a quantum circuit with a given precondition.
    Ex {
        /// The precondition file (HSL format)
        #[arg(value_name = "pre.hsl")]
        pre: String,

        /// The OpenQASM 2.0 circuit file
        #[arg(value_name = "circuit.qasm")]
        circuit: String,
    },

    /// Verify the execution result against a given postcondition.
    Ver {
        /// The precondition file (HSL format)
        #[arg(value_name = "pre.hsl")]
        pre: String,

        /// The OpenQASM 2.0 circuit file
        #[arg(value_name = "circuit.qasm")]
        circuit: String,

        /// The postcondition file (HSL format)
        #[arg(value_name = "post.hsl")]
        post: String,
    },

    /// Check equivalence of two quantum circuits.
    Eq {
        /// The first OpenQASM 2.0 circuit file
        #[arg(value_name = "circuit1.qasm")]
        circuit1: String,

        /// The second OpenQASM 2.0 circuit file
        #[arg(value_name = "circuit2.qasm")]
        circuit2: String,
    },

    /// Print the set of quantum states from an HSL file.
    Print {
        /// The HSL file describing quantum states
        #[arg(value_name = "states.hsl")]
        states: String,
    },
}

fn read_hsl(path: &str) -> Result<Automata<ConcreteSymbol>> {
    let content = fs::read_to_string(path).with_context(|| format!("Failed to read {}", path))?;
    hsl::parse_and_build(&content)
        .map_err(|e| anyhow::anyhow!("HSL parse error in {}: {}", path, e))
}

fn read_qasm(path: &str) -> Result<String> {
    fs::read_to_string(path).with_context(|| format!("Failed to read {}", path))
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let start = Instant::now();

    match cli.command {
        Commands::Ex { pre, circuit } => {
            let mut aut = read_hsl(&pre)?;
            let qasm_src = read_qasm(&circuit)?;

            let gate_count = qasm::execute(&mut aut, &qasm_src)
                .map_err(|e| anyhow::anyhow!("Execution error: {}", e))?;

            println!("OUTPUT:");
            aut.print_aut("  ");
            println!();
            println!("Serialized (Timbuk):");
            print!("{}", timbuk::serialize(&aut));
            println!();

            let elapsed = start.elapsed();
            println!(
                "The quantum program has [{}] qubits and [{}] gates. Executed in [{:.3}s].",
                aut.qubit_num,
                gate_count,
                elapsed.as_secs_f64()
            );
        }

        Commands::Ver {
            pre,
            circuit,
            post,
        } => {
            let mut aut = read_hsl(&pre)?;
            let spec = read_hsl(&post)?;
            let qasm_src = read_qasm(&circuit)?;

            let gate_count = qasm::execute(&mut aut, &qasm_src)
                .map_err(|e| anyhow::anyhow!("Execution error: {}", e))?;

            aut.reduce();
            let mut spec = spec;
            spec.reduce();

            let verify = inclusion::is_included_in(&aut, &spec);
            let elapsed = start.elapsed();

            println!(
                "The quantum program has [{}] qubits and [{}] gates. \
                 The verification process [{}] in [{:.3}s].",
                aut.qubit_num,
                gate_count,
                if verify { "OK" } else { "failed" },
                elapsed.as_secs_f64()
            );
        }

        Commands::Eq { circuit1, circuit2 } => {
            let qasm1 = read_qasm(&circuit1)?;
            let qasm2 = read_qasm(&circuit2)?;

            let n1 = qasm::extract_qubit_count(&qasm1)
                .map_err(|e| anyhow::anyhow!("circuit1: {}", e))?;
            let n2 = qasm::extract_qubit_count(&qasm2)
                .map_err(|e| anyhow::anyhow!("circuit2: {}", e))?;

            if n1 != n2 {
                anyhow::bail!(
                    "Circuits have different qubit counts: {} vs {}",
                    n1,
                    n2
                );
            }

            let mut aut1 = Automata::<ConcreteSymbol>::zero_state(n1);
            let mut aut2 = Automata::<ConcreteSymbol>::zero_state(n2);

            qasm::execute(&mut aut1, &qasm1)
                .map_err(|e| anyhow::anyhow!("circuit1 execution: {}", e))?;
            qasm::execute(&mut aut2, &qasm2)
                .map_err(|e| anyhow::anyhow!("circuit2 execution: {}", e))?;

            aut1.reduce();
            aut2.reduce();

            let result = inclusion::are_equal(&aut1, &aut2);
            let elapsed = start.elapsed();

            println!(
                "The two quantum programs are verified to be [{}] in [{:.3}s].",
                if result { "equal" } else { "unequal" },
                elapsed.as_secs_f64()
            );
        }

        Commands::Print { states } => {
            let aut = read_hsl(&states)?;
            aut.print_aut("");
            println!();
            println!("Serialized (Timbuk):");
            print!("{}", timbuk::serialize(&aut));
        }
    }

    Ok(())
}
