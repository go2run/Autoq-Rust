//! # Timbuk format serialization/deserialization
//!
//! ## C++ correspondence
//! - Serialize: `src/timbuk_serializer.cc`
//! - Deserialize: `src/timbuk_parser-nobison.cc`
//!
//! ## Format
//! ```text
//! Ops
//!
//! Automaton A
//! States q0 q1 q2
//! Final States q0
//! Transitions
//! [1,0,0,0,0](q2, q3) -> q1
//! x1(q0, q1) -> q0
//! ```

use autoq_core::automata::{State, SymbolTag, Tag};
use autoq_core::fivetuple::FiveTuple;
use autoq_core::symbol::ConcreteSymbol;
use autoq_core::Automata;
use num_bigint::BigInt;
use std::collections::BTreeSet;
use std::str::FromStr;

/// Serialize an automaton to Timbuk format string.
pub fn serialize(aut: &Automata<ConcreteSymbol>) -> String {
    let mut result = String::new();

    result.push_str("Ops\n\n");
    result.push_str(&format!(
        "Automaton {}\n",
        if aut.name.is_empty() { "A" } else { &aut.name }
    ));

    // Collect all states
    let mut states: BTreeSet<State> = BTreeSet::new();
    for fs in &aut.final_states {
        states.insert(*fs);
    }
    for (_, map) in &aut.transitions {
        for (parent, children_set) in map {
            states.insert(*parent);
            for children in children_set {
                for &child in children {
                    states.insert(child);
                }
            }
        }
    }

    result.push_str("States");
    for s in &states {
        result.push_str(&format!(" q{}", s));
    }
    result.push('\n');

    result.push_str("Final States");
    for s in &aut.final_states {
        result.push_str(&format!(" q{}", s));
    }
    result.push('\n');

    // Determine if all tags are 1 (single color)
    let all_single_color = aut.transitions.keys().all(|st| st.tag == 1);

    result.push_str("Transitions\n");
    for (st, map) in &aut.transitions {
        let symbol_str = match &st.symbol {
            ConcreteSymbol::Internal(q) => format!("x{}", q),
            ConcreteSymbol::Leaf(ft) => {
                format!("[{},{},{},{},{}]", ft.a, ft.b, ft.c, ft.d, ft.k)
            }
        };

        let tag_str = if !all_single_color && st.tag != 1 {
            let mut bits = Vec::new();
            let mut t = st.tag;
            let mut i = 0u32;
            while t > 0 {
                if t & 1 != 0 {
                    bits.push(i);
                }
                t >>= 1;
                i += 1;
            }
            format!(
                "{{{}}}",
                bits.iter()
                    .map(|b| b.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        } else if !all_single_color && st.tag == 1 {
            "{0}".to_string()
        } else {
            String::new()
        };

        for (parent, children_set) in map {
            for children in children_set {
                if children.is_empty() {
                    result.push_str(&format!("{}{} -> q{}\n", symbol_str, tag_str, parent));
                } else {
                    let children_str = children
                        .iter()
                        .map(|c| format!("q{}", c))
                        .collect::<Vec<_>>()
                        .join(", ");
                    result.push_str(&format!(
                        "{}{}({}) -> q{}\n",
                        symbol_str, tag_str, children_str, parent
                    ));
                }
            }
        }
    }

    result
}

/// Deserialize an automaton from Timbuk format string.
pub fn deserialize(input: &str) -> Result<Automata<ConcreteSymbol>, String> {
    let mut name = String::new();
    let mut final_states: Vec<State> = Vec::new();
    let mut aut = Automata::new(0);
    let mut max_state: State = -1;
    let mut max_qubit: u32 = 0;
    let mut in_transitions = false;

    for line in input.lines() {
        let line = line.trim();

        if line.is_empty() || line == "Ops" {
            continue;
        }

        if let Some(rest) = line.strip_prefix("Automaton ") {
            name = rest.trim().to_string();
            continue;
        }

        if line.starts_with("States") && !line.starts_with("Final States") {
            continue;
        }

        if let Some(rest) = line.strip_prefix("Final States") {
            let rest = rest.trim();
            if !rest.is_empty() {
                for token in rest.split_whitespace() {
                    let s = parse_state(token)?;
                    final_states.push(s);
                    if s > max_state {
                        max_state = s;
                    }
                }
            }
            continue;
        }

        if line == "Transitions" {
            in_transitions = true;
            continue;
        }

        if in_transitions && !line.is_empty() {
            let (tag, rest) = parse_tag_prefix(line)?;
            let (symbol, rest) = parse_symbol(rest)?;
            let (children, rest) = parse_children(rest)?;
            let parent = parse_arrow_target(rest)?;

            if parent > max_state {
                max_state = parent;
            }
            for &c in &children {
                if c > max_state {
                    max_state = c;
                }
            }

            if let ConcreteSymbol::Internal(q) = &symbol {
                if (*q as u32) > max_qubit {
                    max_qubit = *q as u32;
                }
            }

            let st = SymbolTag::new(symbol, tag);
            aut.transitions
                .entry(st)
                .or_default()
                .entry(parent)
                .or_default()
                .insert(children);
        }
    }

    aut.name = name;
    aut.final_states = final_states;
    aut.state_num = max_state + 1;
    aut.qubit_num = max_qubit;

    Ok(aut)
}

fn parse_state(token: &str) -> Result<State, String> {
    let token = token.trim().trim_end_matches(',');
    if let Some(rest) = token.strip_prefix('q') {
        rest.parse::<State>()
            .map_err(|e| format!("Invalid state '{}': {}", token, e))
    } else {
        Err(format!(
            "Expected state starting with 'q', got '{}'",
            token
        ))
    }
}

fn parse_tag_prefix(line: &str) -> Result<(Tag, &str), String> {
    if line.starts_with('{') {
        let end = line
            .find('}')
            .ok_or_else(|| format!("Unclosed '{{' in line: {}", line))?;
        let bits_str = &line[1..end];
        let mut tag: Tag = 0;
        for bit_str in bits_str.split(',') {
            let bit: u32 = bit_str
                .trim()
                .parse()
                .map_err(|e| format!("Invalid tag bit '{}': {}", bit_str, e))?;
            tag |= 1u64 << bit;
        }
        Ok((tag, &line[end + 1..]))
    } else {
        Ok((1, line))
    }
}

fn parse_symbol(s: &str) -> Result<(ConcreteSymbol, &str), String> {
    let s = s.trim();
    if s.starts_with('[') {
        let end = s
            .find(']')
            .ok_or_else(|| format!("Unclosed '[' in: {}", s))?;
        let components_str = &s[1..end];
        let parts: Vec<&str> = components_str.split(',').collect();
        if parts.len() != 5 {
            return Err(format!(
                "FiveTuple needs 5 components, got {}: {}",
                parts.len(),
                components_str
            ));
        }
        let a = BigInt::from_str(parts[0].trim())
            .map_err(|e| format!("Invalid BigInt '{}': {}", parts[0], e))?;
        let b = BigInt::from_str(parts[1].trim())
            .map_err(|e| format!("Invalid BigInt '{}': {}", parts[1], e))?;
        let c = BigInt::from_str(parts[2].trim())
            .map_err(|e| format!("Invalid BigInt '{}': {}", parts[2], e))?;
        let d = BigInt::from_str(parts[3].trim())
            .map_err(|e| format!("Invalid BigInt '{}': {}", parts[3], e))?;
        let k = BigInt::from_str(parts[4].trim())
            .map_err(|e| format!("Invalid BigInt '{}': {}", parts[4], e))?;
        let ft = FiveTuple { a, b, c, d, k };
        Ok((ConcreteSymbol::Leaf(ft), &s[end + 1..]))
    } else if s.starts_with('x') {
        let rest = &s[1..];
        let end = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        let qubit: i64 = rest[..end]
            .parse()
            .map_err(|e| format!("Invalid qubit index in '{}': {}", s, e))?;
        Ok((ConcreteSymbol::Internal(qubit), &rest[end..]))
    } else {
        Err(format!(
            "Expected symbol starting with '[' or 'x', got: {}",
            s
        ))
    }
}

fn parse_children(s: &str) -> Result<(Vec<State>, &str), String> {
    let s = s.trim();
    if s.starts_with('(') {
        let end = s
            .find(')')
            .ok_or_else(|| format!("Unclosed '(' in: {}", s))?;
        let inner = &s[1..end];
        let mut children = Vec::new();
        if !inner.trim().is_empty() {
            for token in inner.split(',') {
                children.push(parse_state(token.trim())?);
            }
        }
        Ok((children, &s[end + 1..]))
    } else {
        Ok((vec![], s))
    }
}

fn parse_arrow_target(s: &str) -> Result<State, String> {
    let s = s.trim();
    let rest = s
        .strip_prefix("->")
        .ok_or_else(|| format!("Expected '->' in: {}", s))?;
    parse_state(rest.trim())
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serialize_zero_state() {
        let aut = Automata::<ConcreteSymbol>::zero_state(2);
        let output = serialize(&aut);
        assert!(output.contains("Ops"));
        assert!(output.contains("Automaton "));
        assert!(output.contains("States"));
        assert!(output.contains("Final States"));
        assert!(output.contains("Transitions"));
        assert!(output.contains("[1,0,0,0,0]"));
        assert!(output.contains("[0,0,0,0,0]"));
        assert!(output.contains("x1"));
        assert!(output.contains("x2"));
    }

    #[test]
    fn serialize_deserialize_roundtrip_zero_state() {
        let original = Automata::<ConcreteSymbol>::zero_state(2);
        let serialized = serialize(&original);
        let deserialized = deserialize(&serialized).expect("deserialization failed");
        assert_eq!(deserialized.final_states.len(), original.final_states.len());
        assert_eq!(deserialized.qubit_num, original.qubit_num);
        assert_eq!(
            deserialized.count_transitions(),
            original.count_transitions()
        );
        let reserialized = serialize(&deserialized);
        assert_eq!(serialized, reserialized);
    }

    #[test]
    fn serialize_uniform_state() {
        let aut = Automata::<ConcreteSymbol>::uniform(2);
        let output = serialize(&aut);
        assert!(output.contains("[1,0,0,0,2]"));
    }

    #[test]
    fn serialize_deserialize_roundtrip_uniform() {
        let original = Automata::<ConcreteSymbol>::uniform(3);
        let serialized = serialize(&original);
        let deserialized = deserialize(&serialized).expect("deserialization failed");
        assert_eq!(deserialized.qubit_num, original.qubit_num);
        assert_eq!(
            deserialized.count_transitions(),
            original.count_transitions()
        );
        let reserialized = serialize(&deserialized);
        assert_eq!(serialized, reserialized);
    }

    #[test]
    fn deserialize_minimal() {
        let input = "\
Ops

Automaton Test
States q0 q1
Final States q0
Transitions
[1,0,0,0,0] -> q1
x1(q1, q1) -> q0
";
        let aut = deserialize(input).expect("deserialization failed");
        assert_eq!(aut.name, "Test");
        assert_eq!(aut.final_states, vec![0]);
        assert_eq!(aut.qubit_num, 1);
        assert_eq!(aut.count_transitions(), 2);
    }

    #[test]
    fn deserialize_with_tags() {
        let input = "\
Ops

Automaton Tagged
States q0 q1
Final States q0
Transitions
{0,1}[1,0,0,0,0] -> q1
{0}x1(q1, q1) -> q0
";
        let aut = deserialize(input).expect("deserialization failed");
        assert_eq!(aut.count_transitions(), 2);
        let has_tag_3 = aut.transitions.keys().any(|st| st.tag == 3);
        assert!(has_tag_3, "expected tag with bits {{0,1}} = 3");
        let has_tag_1 = aut.transitions.keys().any(|st| st.tag == 1);
        assert!(has_tag_1, "expected tag with bit {{0}} = 1");
    }
}
