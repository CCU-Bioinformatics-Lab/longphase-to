#!/usr/bin/env bash
# Check that MethylXGB tuning options are range-checked only when --methyl-xgb
# is given. Every invocation stops during option parsing (no input files), so
# no data is needed; only stderr is inspected.
# Usage: tests/test_phase_option_validation.sh   (after `make`; LONGPHASE overrides the binary)
set -uo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${LONGPHASE:-$root/longphase-to}"
fail=0

run_phase() {
    "$bin" phase "$@" 2>&1 >/dev/null
}

expect() {
    local want="$1" pattern="$2"
    shift 2
    local err
    err="$(run_phase "$@")"
    if grep -qF -- "$pattern" <<<"$err"; then
        [ "$want" = present ] && return 0
    else
        [ "$want" = absent ] && return 0
    fi
    echo "FAIL: expected '$pattern' $want for: phase $*"
    fail=1
}

ignored="ignored because --methyl-xgb is not set"

# MethylXGB off: out-of-range tuning values are not errors, but are reported as ignored.
expect absent "invalid methyl-window"           --methyl-window=0
expect absent "invalid methyl-xgb-snv-threshold" --methyl-xgb-snv-threshold=2
expect absent "invalid methyl-xgb-indel-threshold" --methyl-xgb-indel-threshold=-1
expect absent "invalid methyl thresholds"       --meth-low=0.9 --meth-high=0.1
expect present "$ignored"                       --methyl-window=0
expect absent "invalid methyl-window"           --disable-methyl-xgb --methyl-window=0

# MethylXGB off, no tuning option given: no warning.
expect absent "$ignored"

# Malformed values are still input errors.
expect present "invalid methyl-window"          --methyl-window=abc

# MethylXGB on: out-of-range values are errors and nothing is ignored.
expect present "invalid methyl-window"           --methyl-xgb --methyl-window=0
expect present "invalid methyl-xgb-snv-threshold" --methyl-xgb --methyl-xgb-snv-threshold=2
expect present "invalid methyl-xgb-indel-threshold" --methyl-xgb --methyl-xgb-indel-threshold=-1
expect present "invalid methyl thresholds"       --methyl-xgb --meth-low=0.9 --meth-high=0.1
expect absent "$ignored"                         --methyl-xgb --methyl-window=0

if [ "$fail" -ne 0 ]; then
    exit 1
fi
echo "PASS: phase option validation"
