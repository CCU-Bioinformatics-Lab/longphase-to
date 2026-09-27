#!/usr/bin/env bash
# Check that MethylXGB tuning options are range-checked only when --methyl-xgb
# is given. Uses a header-only VCF, a one-contig reference and an empty BAM, so
# no data is needed. An invocation counts as "accepted" when phase gets past
# option parsing and prints its parameter summary, and as "rejected" when it
# stops at option parsing with a nonzero exit status. The later exit status of
# an accepted run is not checked: these inputs contain no reads.
# Usage: tests/test_phase_option_validation.sh   (after `make`; LONGPHASE overrides the binary)
set -uo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${LONGPHASE:-$root/longphase-to}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n' > "$work/in.vcf"
printf '>chr1\nACGTACGTAC\n' > "$work/ref.fa"
: > "$work/in.bam"

accepted_marker="--- File Parameter ---"
ignored="ignored because --methyl-xgb is not set"
fail=0

# check <accepted|rejected> <pattern> <present|absent> [phase options...]
check() {
    local want="$1" pattern="$2" presence="$3"
    shift 3
    local err status
    err="$(cd "$work" && timeout 60 "$bin" phase -s in.vcf -r ref.fa -b in.bam \
        -c clairs_to_ss -o out "$@" 2>&1 >/dev/null)"
    status=$?
    if [ "$want" = accepted ]; then
        if ! grep -qF -- "$accepted_marker" <<<"$err"; then
            echo "FAIL: expected options accepted for: phase $*"
            fail=1
        fi
    else
        if [ "$status" -eq 0 ] || grep -qF -- "$accepted_marker" <<<"$err"; then
            echo "FAIL: expected options rejected (status $status) for: phase $*"
            fail=1
        fi
    fi
    if grep -qF -- "$pattern" <<<"$err"; then
        [ "$presence" = present ] && return 0
    else
        [ "$presence" = absent ] && return 0
    fi
    echo "FAIL: expected '$pattern' $presence for: phase $*"
    fail=1
}

# MethylXGB off: out-of-range tuning values are accepted and reported as ignored.
check accepted "$ignored" present --methyl-window=0
check accepted "$ignored" present --methyl-xgb-snv-threshold=2
check accepted "$ignored" present --methyl-xgb-indel-threshold=-1
check accepted "$ignored" present --meth-low=0.9 --meth-high=0.1
check accepted "$ignored" present --disable-methyl-xgb --methyl-window=0

# MethylXGB off, no tuning option given: no warning.
check accepted "$ignored" absent

# Malformed values are still input errors, without the ignored warning.
check rejected "invalid methyl-window" present --methyl-window=abc
check rejected "$ignored" absent --methyl-window=abc

# MethylXGB on: out-of-range values are errors.
check rejected "invalid methyl-window" present --methyl-xgb --methyl-window=0
check rejected "invalid methyl-xgb-snv-threshold" present --methyl-xgb --methyl-xgb-snv-threshold=2
check rejected "invalid methyl-xgb-indel-threshold" present --methyl-xgb --methyl-xgb-indel-threshold=-1
check rejected "invalid methyl thresholds" present --methyl-xgb --meth-low=0.9 --meth-high=0.1
check rejected "$ignored" absent --methyl-xgb --methyl-window=0

if [ "$fail" -ne 0 ]; then
    exit 1
fi
echo "PASS: phase option validation"
