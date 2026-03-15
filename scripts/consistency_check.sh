#!/bin/bash
# ──────────────────────────────────────────────────────────────────────────────
# consistency_check.sh — Compare C++ AutoQ and Rust AutoQ-Rust test results
# ──────────────────────────────────────────────────────────────────────────────
#
# This script runs both the C++ AutoQ unit tests and the Rust AutoQ-Rust
# consistency tests, then compares results to verify behavioral equivalence.
#
# Usage:
#   ./scripts/consistency_check.sh [--verbose]
#
# Prerequisites:
#   - C++ AutoQ built at /home/user/AutoQ/build/
#   - Rust AutoQ-Rust workspace at /home/user/Autoq-Rust/
# ──────────────────────────────────────────────────────────────────────────────

set -uo pipefail

VERBOSE="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUST_DIR="$(dirname "$SCRIPT_DIR")"
CPP_DIR="/home/user/AutoQ"
CPP_TEST_BIN="$CPP_DIR/build/unit_tests/explicit_tree_aut_test"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

pass_count=0
fail_count=0
skip_count=0

print_header() {
    echo ""
    echo -e "${BLUE}══════════════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}  AutoQ Consistency Check: C++ vs Rust${NC}"
    echo -e "${BLUE}══════════════════════════════════════════════════════════════${NC}"
    echo ""
}

print_result() {
    local status=$1  # PASS, FAIL, SKIP
    local test_name=$2
    local detail="${3:-}"

    case $status in
        PASS)
            echo -e "  ${GREEN}[PASS]${NC} $test_name"
            ((pass_count++))
            ;;
        FAIL)
            echo -e "  ${RED}[FAIL]${NC} $test_name"
            if [ -n "$detail" ]; then
                echo -e "         ${RED}$detail${NC}"
            fi
            ((fail_count++))
            ;;
        SKIP)
            echo -e "  ${YELLOW}[SKIP]${NC} $test_name"
            if [ -n "$detail" ]; then
                echo -e "         ${YELLOW}$detail${NC}"
            fi
            ((skip_count++))
            ;;
    esac
}

# ── Section 1: Check prerequisites ──────────────────────────────────────────

print_header

echo -e "${BLUE}[1/4] Checking prerequisites...${NC}"

if [ ! -f "$CPP_TEST_BIN" ]; then
    echo -e "  ${RED}C++ test binary not found at $CPP_TEST_BIN${NC}"
    echo "  Run: cd $CPP_DIR && mkdir -p build && cd build && cmake .. && make"
    exit 1
fi
echo "  C++ test binary: OK"

if ! command -v cargo &> /dev/null; then
    echo -e "  ${RED}cargo not found${NC}"
    exit 1
fi
echo "  Rust cargo: OK"

if [ ! -f "$RUST_DIR/Cargo.toml" ]; then
    echo -e "  ${RED}Rust workspace not found at $RUST_DIR${NC}"
    exit 1
fi
echo "  Rust workspace: OK"
echo ""

# ── Section 2: Run C++ gate identity tests ──────────────────────────────────

echo -e "${BLUE}[2/4] Running C++ AutoQ gate identity tests...${NC}"

declare -A CPP_RESULTS

CPP_TESTS=(
    "X_gate_twice_to_identity"
    "Y_gate_twice_to_identity"
    "Z_gate_twice_to_identity"
    "H_gate_twice_to_identity"
    "S_gate_fourth_to_identity"
    "Sdg_gate_equal_to_S_three_times"
    "T_gate_eighth_to_identity"
    "Tdg_gate_equal_to_T_seven_times"
    "swap_gate_simply_exchanges_basis"
    "CX_gate_twice_to_identity"
    "CZ_gate_twice_to_identity"
    "CCX_gate_twice_to_identity"
)

for test in "${CPP_TESTS[@]}"; do
    if $CPP_TEST_BIN --run_test="$test" --log_level=nothing 2>/dev/null; then
        CPP_RESULTS[$test]="PASS"
        print_result "PASS" "C++ $test"
    else
        CPP_RESULTS[$test]="FAIL"
        print_result "FAIL" "C++ $test"
    fi
done

echo ""

# ── Section 3: Run Rust consistency tests ───────────────────────────────────

echo -e "${BLUE}[3/4] Running Rust AutoQ-Rust consistency tests...${NC}"

cd "$RUST_DIR"

# Build first
if ! cargo build -p autoq-core 2>/dev/null; then
    echo -e "  ${RED}Rust build failed!${NC}"
    exit 1
fi

# Run all consistency tests, capture individual results
RUST_OUTPUT=$(cargo test -p autoq-core --test consistency 2>&1)
RUST_EXIT=$?

if [ "$VERBOSE" = "--verbose" ]; then
    echo "$RUST_OUTPUT"
fi

# Parse Rust test results
declare -A RUST_RESULTS

while IFS= read -r line; do
    if [[ "$line" =~ ^test\ ([a-z_]+)\ \.\.\.\ (ok|FAILED) ]]; then
        test_name="${BASH_REMATCH[1]}"
        result="${BASH_REMATCH[2]}"
        if [ "$result" = "ok" ]; then
            RUST_RESULTS[$test_name]="PASS"
            print_result "PASS" "Rust $test_name"
        else
            RUST_RESULTS[$test_name]="FAIL"
            print_result "FAIL" "Rust $test_name"
        fi
    fi
done <<< "$RUST_OUTPUT"

echo ""

# ── Section 4: Cross-comparison ─────────────────────────────────────────────

echo -e "${BLUE}[4/4] Cross-comparison: C++ vs Rust equivalent tests...${NC}"
echo ""

# Map C++ test names to equivalent Rust tests
declare -A CPP_TO_RUST=(
    ["X_gate_twice_to_identity"]="x_gate_twice_identity_zero_state"
    ["Y_gate_twice_to_identity"]="y_gate_twice_identity"
    ["Z_gate_twice_to_identity"]="z_gate_twice_identity"
    ["H_gate_twice_to_identity"]="h_gate_twice_identity"
    ["S_gate_fourth_to_identity"]="s_gate_fourth_identity"
    ["Sdg_gate_equal_to_S_three_times"]="sdg_equals_s_cubed"
    ["T_gate_eighth_to_identity"]="t_gate_eighth_identity"
    ["Tdg_gate_equal_to_T_seven_times"]="tdg_equals_t_seventh"
    ["swap_gate_simply_exchanges_basis"]="swap_gate_identity_on_uniform"
    ["CX_gate_twice_to_identity"]="cx_gate_twice_identity"
    ["CZ_gate_twice_to_identity"]="cz_gate_twice_identity"
    ["CCX_gate_twice_to_identity"]=""  # No Rust equivalent (CCX not implemented)
)

echo "  ┌───────────────────────────────────────────────────────────────┐"
echo "  │ Property              │ C++    │ Rust   │ Match?             │"
echo "  ├───────────────────────────────────────────────────────────────┤"

consistent=0
inconsistent=0
not_comparable=0

for cpp_test in "${CPP_TESTS[@]}"; do
    rust_test="${CPP_TO_RUST[$cpp_test]:-}"
    cpp_result="${CPP_RESULTS[$cpp_test]:-N/A}"

    if [ -z "$rust_test" ]; then
        printf "  │ %-21s │ %-6s │ %-6s │ %-18s │\n" \
            "$cpp_test" "$cpp_result" "N/A" "Not implemented"
        ((not_comparable++))
        continue
    fi

    rust_result="${RUST_RESULTS[$rust_test]:-N/A}"

    if [ "$cpp_result" = "$rust_result" ]; then
        match="YES"
        ((consistent++))
    else
        match="DIFFERS"
        ((inconsistent++))
    fi

    # Shorten names for table
    short_name=$(echo "$cpp_test" | sed 's/_to_identity//' | sed 's/_gate//' | head -c 21)
    printf "  │ %-21s │ %-6s │ %-6s │ %-18s │\n" \
        "$short_name" "$cpp_result" "$rust_result" "$match"
done

echo "  └───────────────────────────────────────────────────────────────┘"
echo ""

# ── Summary ─────────────────────────────────────────────────────────────────

echo -e "${BLUE}══════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}  Summary${NC}"
echo -e "${BLUE}══════════════════════════════════════════════════════════════${NC}"
echo ""

# Rust-only tests summary
rust_total=${#RUST_RESULTS[@]}
rust_pass=0
rust_fail=0
for test in "${!RUST_RESULTS[@]}"; do
    if [ "${RUST_RESULTS[$test]}" = "PASS" ]; then
        ((rust_pass++))
    else
        ((rust_fail++))
    fi
done

echo "  C++ AutoQ tests:     ${#CPP_RESULTS[@]} total"
echo "  Rust AutoQ-Rust:     $rust_total total ($rust_pass pass, $rust_fail fail)"
echo ""
echo "  Cross-comparison:"
echo -e "    Consistent:      ${GREEN}$consistent${NC}"
echo -e "    Inconsistent:    ${RED}$inconsistent${NC}"
echo -e "    Not comparable:  ${YELLOW}$not_comparable${NC} (not yet implemented in Rust)"
echo ""

# Rust-only tests (no C++ equivalent)
echo "  Rust-only tests (no C++ equivalent):"
for test in "${!RUST_RESULTS[@]}"; do
    is_mapped=false
    for rust_test in "${CPP_TO_RUST[@]}"; do
        if [ "$test" = "$rust_test" ]; then
            is_mapped=true
            break
        fi
    done
    if ! $is_mapped; then
        result="${RUST_RESULTS[$test]}"
        if [ "$result" = "PASS" ]; then
            echo -e "    ${GREEN}[PASS]${NC} $test"
        else
            echo -e "    ${RED}[FAIL]${NC} $test"
        fi
    fi
done

echo ""

# Known differences
echo "  Known differences (expected):"
echo "    - C++ X/CX/CCX/Swap tests include 'random' automaton which can fail"
echo "      due to nondeterminism in the inclusion check."
echo "    - Rust uses 3 qubits (vs C++ 14) for faster testing."
echo "    - Rust HXH=Z/HZH=X on zero_state deferred due to inclusion limitation."
echo "    - CCX (Toffoli) not yet implemented in Rust."
echo ""

if [ $inconsistent -eq 0 ]; then
    echo -e "  ${GREEN}Result: All comparable tests are CONSISTENT between C++ and Rust.${NC}"
    exit 0
else
    echo -e "  ${RED}Result: $inconsistent test(s) show INCONSISTENT behavior.${NC}"
    exit 1
fi
