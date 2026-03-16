//! # HSL parser — Extended Dirac notation
//!
//! ## C++ correspondence
//! Replaces `src/ExtendedDirac/` (ANTLR4-generated parser) with pest PEG.
//!
//! ## Supported formats
//! - Single-state: `{c1 |0101>}`
//! - Multi-state sets: `{|00>, |01>, 1/sqrt2 |00> + 1/sqrt2 |01>}`
//! - Rational amplitudes: `{75555... / (sqrt2^152) |0000001>}`
//!
//! ## Not yet supported
//! - Symbolic kets with variable bits (`|ss0000001>`)
//! - Sum quantifiers (`∑ |i|=8, i≠s`)
//! - Constrained postconditions (`:` syntax)
//! These require integration with Z3 and a more complex grammar.

use autoq_core::automata::{State, Tag};
use autoq_core::fivetuple::FiveTuple;
use autoq_core::symbol::ConcreteSymbol;
use autoq_core::Automata;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use pest::Parser;
use pest_derive::Parser;
use std::collections::HashMap;

// ─── Pest parser derivation ─────────────────────────────────────────────────

#[derive(Parser)]
#[grammar = "grammar.pest"]
struct HslParser;

// ─── Intermediate AST types ─────────────────────────────────────────────────

/// A parsed amplitude value (before constant substitution).
#[derive(Debug, Clone)]
pub enum AmpExpr {
    /// A named constant reference (e.g. `c1`, `aH`).
    Const(String),
    /// `1/sqrt2`
    InvSqrt2,
    /// `1/2`
    Half,
    /// A rational `num / (sqrt2^exp)`.
    RationalSqrt2 { num: BigInt, sqrt2_exp: u32 },
    /// A rational `num / den`.
    Rational { num: BigInt, den: BigInt },
    /// A bare integer.
    Integer(BigInt),
}

impl AmpExpr {
    /// Evaluate to a concrete FiveTuple, substituting constants from `env`.
    pub fn eval(&self, env: &HashMap<String, FiveTuple>) -> Result<FiveTuple, String> {
        match self {
            AmpExpr::Const(name) => env
                .get(name)
                .cloned()
                .ok_or_else(|| format!("Undefined constant '{}'", name)),
            AmpExpr::InvSqrt2 => Ok(FiveTuple::inv_sqrt2()),
            AmpExpr::Half => Ok(FiveTuple {
                a: BigInt::one(),
                b: BigInt::zero(),
                c: BigInt::zero(),
                d: BigInt::zero(),
                k: BigInt::from(2u32), // 1/√2^2 = 1/2
            }),
            AmpExpr::RationalSqrt2 { num, sqrt2_exp } => Ok(FiveTuple {
                a: num.clone(),
                b: BigInt::zero(),
                c: BigInt::zero(),
                d: BigInt::zero(),
                k: BigInt::from(*sqrt2_exp),
            }),
            AmpExpr::Rational { num, den } => {
                let mut d = den.clone();
                let mut k: u32 = 0;
                let two = BigInt::from(2u32);
                while &d % &two == BigInt::zero() {
                    d /= &two;
                    k += 2;
                }
                if d != BigInt::one() {
                    return Err(format!(
                        "Denominator {} is not a power of 2 (cannot represent in FiveTuple)",
                        den
                    ));
                }
                Ok(FiveTuple {
                    a: num.clone(),
                    b: BigInt::zero(),
                    c: BigInt::zero(),
                    d: BigInt::zero(),
                    k: BigInt::from(k),
                })
            }
            AmpExpr::Integer(n) => Ok(FiveTuple::from_int(n.clone())),
        }
    }
}

/// A single `amplitude * |ket⟩` term.
#[derive(Debug, Clone)]
pub struct Term {
    pub sign: i32, // +1 or -1
    pub amplitude: AmpExpr,
    pub ket: Vec<u8>, // bit string: 0s and 1s
}

/// A quantum state: a linear combination of terms.
#[derive(Debug, Clone)]
pub struct StateExpr {
    pub terms: Vec<Term>,
}

/// The full parsed HSL file.
#[derive(Debug)]
pub struct HslFile {
    pub constants: HashMap<String, AmpExpr>,
    pub states: Vec<StateExpr>,
}

// ─── Parsing logic ──────────────────────────────────────────────────────────

pub fn parse_hsl(input: &str) -> Result<HslFile, String> {
    let pairs = HslParser::parse(Rule::file, input).map_err(|e| format!("Parse error: {}", e))?;

    let mut constants: HashMap<String, AmpExpr> = HashMap::new();
    let mut states: Vec<StateExpr> = Vec::new();

    for pair in pairs {
        match pair.as_rule() {
            Rule::file => {
                for inner in pair.into_inner() {
                    match inner.as_rule() {
                        Rule::constants_section => {
                            for def in inner.into_inner() {
                                if def.as_rule() == Rule::const_def {
                                    let (name, expr) = parse_const_def(def)?;
                                    constants.insert(name, expr);
                                }
                            }
                        }
                        Rule::dirac_section => {
                            for dirac_child in inner.into_inner() {
                                match dirac_child.as_rule() {
                                    Rule::dirac_set | Rule::dirac_bare => {
                                        for state_pair in dirac_child.into_inner() {
                                            if state_pair.as_rule() == Rule::state_expr {
                                                states.push(parse_state_expr(state_pair)?);
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }

    Ok(HslFile { constants, states })
}

fn parse_const_def(pair: pest::iterators::Pair<Rule>) -> Result<(String, AmpExpr), String> {
    let mut inner = pair.into_inner();
    let name = inner.next().ok_or("missing name")?.as_str().to_string();
    let amp_pair = inner.next().ok_or("missing value")?;
    let expr = parse_amp_value(amp_pair)?;
    Ok((name, expr))
}

fn parse_amp_value(pair: pest::iterators::Pair<Rule>) -> Result<AmpExpr, String> {
    match pair.as_rule() {
        Rule::amp_value => {
            let inner = pair.into_inner().next().ok_or("empty amp_value")?;
            parse_amp_value(inner)
        }
        Rule::rational_amp => {
            let mut parts = pair.into_inner();
            let num_str = parts.next().ok_or("missing numerator")?.as_str();
            let num: BigInt = num_str.parse().map_err(|e| format!("Bad int: {}", e))?;
            let denom_pair = parts.next().ok_or("missing denominator")?;
            match denom_pair.as_rule() {
                Rule::sqrt2_power => {
                    let exp_str = denom_pair
                        .into_inner()
                        .next()
                        .ok_or("missing exp")?
                        .as_str();
                    let exp: u32 = exp_str.parse().map_err(|e| format!("Bad exp: {}", e))?;
                    Ok(AmpExpr::RationalSqrt2 {
                        num,
                        sqrt2_exp: exp,
                    })
                }
                Rule::uint => {
                    let den: BigInt = denom_pair.as_str().parse().map_err(|e| format!("{}", e))?;
                    Ok(AmpExpr::Rational { num, den })
                }
                _ => Err(format!(
                    "Unexpected denominator rule: {:?}",
                    denom_pair.as_rule()
                )),
            }
        }
        Rule::identifier => Ok(AmpExpr::Const(pair.as_str().to_string())),
        Rule::inv_sqrt2_lit => Ok(AmpExpr::InvSqrt2),
        Rule::half_lit => Ok(AmpExpr::Half),
        Rule::uint => {
            let n: BigInt = pair.as_str().parse().map_err(|e| format!("{}", e))?;
            Ok(AmpExpr::Integer(n))
        }
        _ => Err(format!(
            "Unknown amp rule: {:?} text={}",
            pair.as_rule(),
            pair.as_str()
        )),
    }
}

fn parse_state_expr(pair: pest::iterators::Pair<Rule>) -> Result<StateExpr, String> {
    let mut terms = Vec::new();
    let mut current_sign = 1i32;

    for child in pair.into_inner() {
        match child.as_rule() {
            Rule::sign => {
                current_sign = if child.as_str().trim() == "-" { -1 } else { 1 };
            }
            Rule::term => {
                let term = parse_term(child, current_sign)?;
                terms.push(term);
                current_sign = 1;
            }
            _ => {}
        }
    }

    Ok(StateExpr { terms })
}

fn parse_term(pair: pest::iterators::Pair<Rule>, sign: i32) -> Result<Term, String> {
    let mut amplitude: Option<AmpExpr> = None;
    let mut ket: Vec<u8> = Vec::new();

    for child in pair.into_inner() {
        match child.as_rule() {
            Rule::amp_value => {
                amplitude = Some(parse_amp_value(child)?);
            }
            Rule::ket => {
                for bit_pair in child.into_inner() {
                    ket.push(bit_pair.as_str().parse::<u8>().unwrap());
                }
            }
            _ => {}
        }
    }

    let amplitude = amplitude.unwrap_or(AmpExpr::Integer(BigInt::one()));
    Ok(Term {
        sign,
        amplitude,
        ket,
    })
}

// ─── Build Automata from parsed HSL ─────────────────────────────────────────

/// Convert a parsed HSL file into a concrete `Automata<ConcreteSymbol>`.
///
/// For multi-state sets (used in nondeterministic specs), returns an automaton
/// whose language is the union (via additional colors) of all states.
///
/// ## Multi-state encoding
/// The C++ codebase encodes each state in the set with a unique "color" (Tag bit).
/// Transitions belonging to state `i` carry tag `1 << i`.  The automaton accepts
/// a tree if all transitions along the path share a common color bit.
pub fn build_automata(hsl: &HslFile) -> Result<Automata<ConcreteSymbol>, String> {
    if hsl.states.is_empty() {
        return Err("No quantum states in HSL file".into());
    }

    // Evaluate constants to FiveTuples
    let mut const_env: HashMap<String, FiveTuple> = HashMap::new();
    for (name, expr) in &hsl.constants {
        const_env.insert(name.clone(), expr.eval(&const_env)?);
    }

    // Determine qubit count from first state
    let n_qubits = hsl.states[0]
        .terms
        .first()
        .map(|t| t.ket.len() as u32)
        .unwrap_or(0);
    if n_qubits == 0 {
        return Err("Empty ket".into());
    }

    let mut aut = Automata::new(n_qubits);
    aut.name = "from_hsl".to_string();

    // Encode each state in the set with its own color bit
    for (state_idx, state_expr) in hsl.states.iter().enumerate() {
        let color: Tag = 1u64 << state_idx;
        build_state_into_automata(state_expr, &const_env, &mut aut, color, n_qubits)?;
    }

    Ok(aut)
}

/// Encode a single quantum state (linear combination of |ket⟩ terms) into the
/// automaton, using the given `color` tag for all its transitions.
fn build_state_into_automata(
    state: &StateExpr,
    env: &HashMap<String, FiveTuple>,
    aut: &mut Automata<ConcreteSymbol>,
    color: Tag,
    n_qubits: u32,
) -> Result<(), String> {
    let num_basis = 1usize << n_qubits;
    let mut amplitudes: HashMap<Vec<u8>, FiveTuple> = HashMap::new();

    for term in &state.terms {
        if term.ket.len() != n_qubits as usize {
            return Err(format!(
                "Ket length {} doesn't match qubit count {}",
                term.ket.len(),
                n_qubits
            ));
        }
        let mut amp = term.amplitude.eval(env)?;
        if term.sign == -1 {
            amp = -amp;
        }
        let entry = amplitudes
            .entry(term.ket.clone())
            .or_insert_with(FiveTuple::zero);
        *entry = entry.clone() + amp;
    }

    // Create leaf states, deduplicating equal amplitudes
    let mut amp_to_state: HashMap<FiveTuple, State> = HashMap::new();

    let get_or_create_leaf =
        |aut: &mut Automata<ConcreteSymbol>,
         amp_to_state: &mut HashMap<FiveTuple, State>,
         amp: FiveTuple|
         -> State {
            if let Some(&s) = amp_to_state.get(&amp) {
                return s;
            }
            let s = aut.new_state();
            aut.add_transition(ConcreteSymbol::Leaf(amp.clone()), color, s, vec![]);
            amp_to_state.insert(amp, s);
            s
        };

    // Build leaf states for all basis kets
    let mut prefix_to_state: HashMap<Vec<u8>, State> = HashMap::new();

    for ket_bits in 0..num_basis {
        let bits: Vec<u8> = (0..n_qubits)
            .map(|i| ((ket_bits >> i) & 1) as u8)
            .rev()
            .collect();
        let amp = amplitudes.get(&bits).cloned().unwrap_or_else(FiveTuple::zero);
        let leaf_state = get_or_create_leaf(aut, &mut amp_to_state, amp);
        prefix_to_state.insert(bits, leaf_state);
    }

    // Build internal nodes bottom-up
    for qubit in (1..=n_qubits).rev() {
        let mut next_prefix_to_state: HashMap<Vec<u8>, State> = HashMap::new();
        let prefix_len = (qubit - 1) as usize;
        let bit_idx = (qubit - 1) as usize;

        let mut seen: HashMap<Vec<u8>, (State, State)> = HashMap::new();
        for (bits, &child_state) in &prefix_to_state {
            let prefix = bits[..prefix_len].to_vec();
            let bit = bits[bit_idx];
            let entry = seen
                .entry(prefix.clone())
                .or_insert((State::MAX, State::MAX));
            if bit == 0 {
                entry.0 = child_state;
            } else {
                entry.1 = child_state;
            }
        }

        for (prefix, (c0, c1)) in seen {
            let s = aut.new_state();
            aut.add_transition(
                ConcreteSymbol::Internal(qubit as i64),
                color,
                s,
                vec![c0, c1],
            );
            next_prefix_to_state.insert(prefix, s);
        }
        prefix_to_state = next_prefix_to_state;
    }

    // The root state (empty prefix) becomes a final state
    if let Some(&root) = prefix_to_state.get(&vec![]) {
        if !aut.final_states.contains(&root) {
            aut.final_states.push(root);
        }
    }

    Ok(())
}

/// Parse an HSL file and build the automaton in one step.
pub fn parse_and_build(input: &str) -> Result<Automata<ConcreteSymbol>, String> {
    let hsl = parse_hsl(input)?;
    build_automata(&hsl)
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const HSL_SINGLE_STATE: &str = r#"Constants
c1 := 1
Extended Dirac
{c1 |0000>}
"#;

    const HSL_SUPERPOSITION: &str = r#"Extended Dirac
{1/sqrt2 |00> + 1/sqrt2 |11>}
"#;

    const HSL_MULTI_STATE_SET: &str = r#"Extended Dirac
{|00>, |01>, 1/sqrt2 |00> + 1/sqrt2 |01>}
"#;

    const HSL_RATIONAL_AMP: &str = r#"Constants
aH := 75555863006653472909761 / (sqrt2 ^ 152)
Extended Dirac
{aH |10000000>}
"#;

    #[test]
    fn parse_single_state() {
        let hsl = parse_hsl(HSL_SINGLE_STATE).expect("parse failed");
        assert_eq!(hsl.constants.len(), 1);
        assert_eq!(hsl.states.len(), 1);
        assert_eq!(hsl.states[0].terms.len(), 1);
        let term = &hsl.states[0].terms[0];
        assert_eq!(term.ket, vec![0, 0, 0, 0]);
    }

    #[test]
    fn parse_superposition() {
        let hsl = parse_hsl(HSL_SUPERPOSITION).expect("parse failed");
        assert_eq!(hsl.states.len(), 1);
        assert_eq!(hsl.states[0].terms.len(), 2);
    }

    #[test]
    fn parse_multi_state_set() {
        let hsl = parse_hsl(HSL_MULTI_STATE_SET).expect("parse failed");
        assert_eq!(hsl.states.len(), 3);
    }

    #[test]
    fn parse_rational_amplitude() {
        let hsl = parse_hsl(HSL_RATIONAL_AMP).expect("parse failed");
        assert_eq!(hsl.constants.len(), 1);
        let mut env = HashMap::new();
        for (name, expr) in &hsl.constants {
            let ft = expr.eval(&env).expect("eval failed");
            env.insert(name.clone(), ft);
        }
    }

    #[test]
    fn build_single_state_automata() {
        let hsl = parse_hsl(HSL_SINGLE_STATE).expect("parse failed");
        let aut = build_automata(&hsl).expect("build failed");
        assert_eq!(aut.qubit_num, 4);
        assert!(!aut.final_states.is_empty());
    }

    #[test]
    fn build_superposition_automata() {
        let hsl = parse_hsl(HSL_SUPERPOSITION).expect("parse failed");
        let aut = build_automata(&hsl).expect("build failed");
        assert_eq!(aut.qubit_num, 2);
    }

    #[test]
    fn build_multi_state_set() {
        let hsl = parse_hsl(HSL_MULTI_STATE_SET).expect("parse failed");
        let aut = build_automata(&hsl).expect("build failed");
        assert_eq!(aut.qubit_num, 2);
        let tags: std::collections::BTreeSet<u64> =
            aut.transitions.keys().map(|st| st.tag).collect();
        assert!(tags.len() >= 3, "Expected ≥3 colors, got {:?}", tags);
    }

    #[test]
    fn real_benchmark_bv_pre() {
        let hsl_str = "Constants\nc1 := 1\nExtended Dirac\n{c1 |00000000>}\n";
        let hsl = parse_hsl(hsl_str).expect("parse failed");
        let aut = build_automata(&hsl).expect("build failed");
        assert_eq!(aut.qubit_num, 8);
    }

    #[test]
    fn parse_fixture_groverfor_pre() {
        let input = std::fs::read_to_string(
            concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/testcase/GroverFor/pre.hsl")
        ).expect("fixture not found");
        let aut = parse_and_build(&input).expect("parse+build failed");
        assert_eq!(aut.qubit_num, 5);
    }
}
