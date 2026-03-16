//! # QASM parser — OpenQASM 2.0 line-by-line gate executor
//!
//! ## C++ correspondence
//! Port of `src/execute.cc` `single_gate_execute()` — regex-based line parsing.
//!
//! Parses OpenQASM 2.0 files line by line and applies gates to an automaton.
//! Supports: x, y, z, h, s, sdg, t, tdg, rx, rz, ry(pi/2), cx, cz, ccx, swap.

use autoq_core::symbol::SymbolTrait;
use autoq_core::Automata;

/// Extract the qubit register size from a QASM file.
pub fn extract_qubit_count(input: &str) -> Result<u32, String> {
    for line in input.lines() {
        let line = line.trim();
        if line.starts_with("qreg ") || line.starts_with("qubit") {
            if let Some(n) = extract_first_number(line) {
                return Ok(n as u32);
            }
        }
    }
    Err("No qreg declaration found in QASM file".into())
}

/// Execute a QASM circuit on an automaton, applying gates line by line.
///
/// This is the Rust equivalent of C++ `Automata::execute()` + `single_gate_execute()`.
/// The automaton is modified in place.
pub fn execute<S: SymbolTrait>(aut: &mut Automata<S>, qasm: &str) -> Result<u32, String> {
    let mut gate_count: u32 = 0;
    let mut in_gate_def = false;
    let mut in_loop = false;
    let mut loop_body: Vec<String> = Vec::new();
    let mut loop_count = 0u32;

    for line in qasm.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        // Skip gate definitions
        if in_gate_def {
            if line.contains('}') {
                in_gate_def = false;
            }
            continue;
        }

        // Skip headers and comments
        if line.starts_with("OPENQASM")
            || line.starts_with("include ")
            || line.starts_with("//")
            || line.starts_with("/*")
            || line.starts_with("bit")
        {
            continue;
        }

        if line.starts_with("qreg ") || line.starts_with("qubit") {
            // Validate qubit count matches
            if let Some(n) = extract_first_number(line) {
                if n as u32 != aut.qubit_num {
                    return Err(format!(
                        "QASM declares {} qubits but automaton has {}",
                        n, aut.qubit_num
                    ));
                }
            }
            continue;
        }

        if line.starts_with("gate ") {
            if !line.contains('}') {
                in_gate_def = true;
            }
            continue;
        }

        // For loop handling: `for int i in [start:end] {`
        if line.starts_with("for ") {
            if in_loop {
                return Err("Nested loops are not supported".into());
            }
            in_loop = true;
            loop_body.clear();
            // Parse loop bounds: for int i in [start:end] {
            let (start, end) = parse_loop_bounds(line)?;
            loop_count = (end - start + 1) as u32;
            continue;
        }

        if line.starts_with('{') && in_loop {
            continue;
        }

        if line.starts_with('}') {
            if in_loop {
                // Execute loop body `loop_count` times
                for _ in 0..loop_count {
                    for loop_line in &loop_body.clone() {
                        gate_count += execute_single_gate(aut, loop_line)?;
                    }
                }
                in_loop = false;
                loop_body.clear();
            }
            continue;
        }

        if in_loop {
            loop_body.push(line.to_string());
            continue;
        }

        // Control flow / debug commands
        if line.starts_with("PRINT_STATS") || line.starts_with("PRINT_AUT") {
            continue;
        }
        if line.starts_with("STOP") {
            break;
        }

        gate_count += execute_single_gate(aut, line)?;
    }

    Ok(gate_count)
}

/// Execute a single gate line on the automaton.
/// Returns 1 if a gate was applied, 0 otherwise.
fn execute_single_gate<S: SymbolTrait>(aut: &mut Automata<S>, line: &str) -> Result<u32, String> {
    // Note: C++ uses 1-based qubit indexing internally (qubit 1..n),
    // while QASM uses 0-based (qubits[0]..qubits[n-1]).
    // We add 1 to convert.

    if line.starts_with("x ") {
        let q = extract_first_number(line).ok_or("x: missing qubit")? as u32;
        aut.x(1 + q);
    } else if line.starts_with("y ") {
        let q = extract_first_number(line).ok_or("y: missing qubit")? as u32;
        aut.y(1 + q);
    } else if line.starts_with("z ") {
        let q = extract_first_number(line).ok_or("z: missing qubit")? as u32;
        aut.z(1 + q);
    } else if line.starts_with("h ") {
        let q = extract_first_number(line).ok_or("h: missing qubit")? as u32;
        aut.h(1 + q);
    } else if line.starts_with("s ") {
        let q = extract_first_number(line).ok_or("s: missing qubit")? as u32;
        aut.s_gate(1 + q);
    } else if line.starts_with("sdg ") {
        let q = extract_first_number(line).ok_or("sdg: missing qubit")? as u32;
        aut.sdg(1 + q);
    } else if line.starts_with("t ") {
        let q = extract_first_number(line).ok_or("t: missing qubit")? as u32;
        aut.t_gate(1 + q);
    } else if line.starts_with("tdg ") {
        let q = extract_first_number(line).ok_or("tdg: missing qubit")? as u32;
        aut.tdg(1 + q);
    } else if line.starts_with("rx(") {
        let (theta_num, theta_den, qubit) = parse_rotation_gate(line, "rx")?;
        aut.rx(theta_num, theta_den, 1 + qubit);
    } else if line.starts_with("rz(") {
        let (theta_num, theta_den, qubit) = parse_rotation_gate(line, "rz")?;
        aut.rz(theta_num, theta_den, 1 + qubit);
    } else if line.starts_with("ry(pi/2)") || line.starts_with("ry(pi / 2)") {
        let nums = extract_all_numbers(line);
        if nums.len() < 2 {
            return Err("ry: missing qubit".into());
        }
        // nums[0] is 2 from "pi/2", nums[1] is the qubit
        aut.ry(1 + nums[nums.len() - 1] as u32);
    } else if line.starts_with("cx ") || line.starts_with("CX ") {
        let nums = extract_all_numbers(line);
        if nums.len() < 2 {
            return Err("cx: need 2 qubits".into());
        }
        aut.cx(1 + nums[0] as u32, 1 + nums[1] as u32);
    } else if line.starts_with("cz ") {
        let nums = extract_all_numbers(line);
        if nums.len() < 2 {
            return Err("cz: need 2 qubits".into());
        }
        aut.cz(1 + nums[0] as u32, 1 + nums[1] as u32);
    } else if line.starts_with("ccx ") {
        let nums = extract_all_numbers(line);
        if nums.len() < 3 {
            return Err("ccx: need 3 qubits".into());
        }
        aut.ccx(1 + nums[0] as u32, 1 + nums[1] as u32, 1 + nums[2] as u32);
    } else if line.starts_with("swap ") {
        let nums = extract_all_numbers(line);
        if nums.len() < 2 {
            return Err("swap: need 2 qubits".into());
        }
        aut.swap(1 + nums[0] as u32, 1 + nums[1] as u32);
    } else if !line.is_empty() {
        return Err(format!("unsupported gate: {}", line));
    } else {
        return Ok(0);
    }

    // Reduce after each gate to keep automaton size manageable
    aut.reduce();

    Ok(1)
}

/// Extract the first decimal number from a string.
fn extract_first_number(s: &str) -> Option<i64> {
    let mut start = None;
    for (i, c) in s.char_indices() {
        if c.is_ascii_digit() {
            if start.is_none() {
                start = Some(i);
            }
        } else if start.is_some() {
            return s[start.unwrap()..i].parse().ok();
        }
    }
    if let Some(st) = start {
        return s[st..].parse().ok();
    }
    None
}

/// Extract all decimal numbers from a string.
fn extract_all_numbers(s: &str) -> Vec<i64> {
    let mut nums = Vec::new();
    let mut start = None;
    for (i, c) in s.char_indices() {
        if c.is_ascii_digit() {
            if start.is_none() {
                start = Some(i);
            }
        } else if let Some(st) = start {
            if let Ok(n) = s[st..i].parse() {
                nums.push(n);
            }
            start = None;
        }
    }
    if let Some(st) = start {
        if let Ok(n) = s[st..].parse() {
            nums.push(n);
        }
    }
    nums
}

/// Parse loop bounds from `for int i in [start:end] {` or `for i in [start:end]`.
fn parse_loop_bounds(line: &str) -> Result<(i64, i64), String> {
    let bracket_start = line.find('[').ok_or("loop: missing '['")?;
    let bracket_end = line.find(']').ok_or("loop: missing ']'")?;
    let inner = &line[bracket_start + 1..bracket_end];
    let parts: Vec<&str> = inner.split(':').collect();
    if parts.len() != 2 {
        return Err(format!("loop: expected [start:end], got [{}]", inner));
    }
    let start: i64 = parts[0].trim().parse().map_err(|e| format!("loop start: {}", e))?;
    let end: i64 = parts[1].trim().parse().map_err(|e| format!("loop end: {}", e))?;
    Ok((start, end))
}

/// Parse a rotation gate line like `rx(pi/4) qubits[2];` or `rz(3*pi/8) qubits[0];`
///
/// Returns (theta_num, theta_den, qubit) where the angle is θ = num/den * π.
/// The C++ code replaces "pi" with "(1/2)" then evaluates, giving θ/2 as a rational.
/// We instead parse the angle directly: `rx(pi/N)` → theta = 1/N of a full turn,
/// meaning the rotation angle is π/N. For the FiveTuple rotation methods,
/// θ_num/θ_den represents the angle as a fraction of 2π.
/// So `rx(pi/4)` → rotation by π/4 → θ_num=1, θ_den=8 (since π/4 = 2π/8).
fn parse_rotation_gate(line: &str, gate: &str) -> Result<(i64, i64, u32), String> {
    // Extract angle string between parentheses
    let paren_start = line.find('(').ok_or(format!("{}: missing '('", gate))?;
    let paren_end = line.find(')').ok_or(format!("{}: missing ')'", gate))?;
    let angle_str = line[paren_start + 1..paren_end].trim();

    // Extract qubit from the rest
    let rest = &line[paren_end + 1..];
    let qubit = extract_first_number(rest).ok_or(format!("{}: missing qubit", gate))? as u32;

    // Parse angle: support "pi/N", "N*pi/M", "pi", "0"
    let (theta_num, theta_den) = parse_angle(angle_str)?;

    Ok((theta_num, theta_den, qubit))
}

/// Parse an angle expression to (numerator, denominator) where angle = num/den * π.
/// The returned values are suitable for FiveTuple rotation methods which expect
/// the angle as a fraction where counterclockwise(num, den) rotates by num/den of 2π.
/// So π/4 becomes (1, 8) since π/4 = (1/8) * 2π.
fn parse_angle(s: &str) -> Result<(i64, i64), String> {
    let s = s.trim();

    if s == "0" {
        return Ok((0, 1));
    }

    if s == "pi" {
        // π = (1/2) * 2π
        return Ok((1, 2));
    }

    // "pi/N" → π/N = (1/(2N)) * 2π → (1, 2N)
    if let Some(rest) = s.strip_prefix("pi/") {
        let den: i64 = rest.trim().parse().map_err(|e| format!("angle: {}", e))?;
        return Ok((1, 2 * den));
    }

    // "pi / N"
    if s.starts_with("pi") {
        let rest = s.strip_prefix("pi").unwrap().trim();
        if let Some(rest) = rest.strip_prefix('/') {
            let den: i64 = rest.trim().parse().map_err(|e| format!("angle: {}", e))?;
            return Ok((1, 2 * den));
        }
    }

    // "N*pi/M" → N*π/M = (N/(2M)) * 2π
    if s.contains("pi") {
        let parts: Vec<&str> = s.split("pi").collect();
        if parts.len() == 2 {
            let num_part = parts[0].trim().trim_end_matches('*').trim();
            let den_part = parts[1].trim().trim_start_matches('/').trim();
            let num: i64 = if num_part.is_empty() {
                1
            } else {
                num_part.parse().map_err(|e| format!("angle num: {}", e))?
            };
            let den: i64 = if den_part.is_empty() {
                1
            } else {
                den_part.parse().map_err(|e| format!("angle den: {}", e))?
            };
            return Ok((num, 2 * den));
        }
    }

    Err(format!("Cannot parse angle: '{}'", s))
}

#[cfg(test)]
mod tests {
    use super::*;
    use autoq_core::symbol::ConcreteSymbol;

    #[test]
    fn extract_qubit_count_basic() {
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg qubits[3];\n";
        assert_eq!(extract_qubit_count(qasm).unwrap(), 3);
    }

    #[test]
    fn parse_angle_tests() {
        assert_eq!(parse_angle("0").unwrap(), (0, 1));
        assert_eq!(parse_angle("pi").unwrap(), (1, 2));
        assert_eq!(parse_angle("pi/4").unwrap(), (1, 8));
        assert_eq!(parse_angle("pi/2").unwrap(), (1, 4));
        assert_eq!(parse_angle("3*pi/4").unwrap(), (3, 8));
    }

    #[test]
    fn execute_empty_circuit() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(2);
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\n\n";
        let count = execute(&mut aut, qasm).unwrap();
        assert_eq!(count, 0);
    }

    #[test]
    fn execute_x_gate() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nx q[0];\n";
        let count = execute(&mut aut, qasm).unwrap();
        assert_eq!(count, 1);
    }

    #[test]
    fn execute_h_gate() {
        let mut aut = Automata::<ConcreteSymbol>::zero_state(1);
        let qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nh q[0];\n";
        let count = execute(&mut aut, qasm).unwrap();
        assert_eq!(count, 1);
    }

    #[test]
    fn extract_numbers() {
        assert_eq!(extract_all_numbers("cx qubits[0], qubits[1];"), vec![0, 1]);
        assert_eq!(
            extract_all_numbers("ccx q[0], q[1], q[2];"),
            vec![0, 1, 2]
        );
    }

    #[test]
    fn loop_bounds() {
        assert_eq!(parse_loop_bounds("for int i in [1:2] {").unwrap(), (1, 2));
        assert_eq!(parse_loop_bounds("for i in [0:3]").unwrap(), (0, 3));
    }
}
