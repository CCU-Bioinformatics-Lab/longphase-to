#!/usr/bin/env bash
# Check the haplotag --log row of a read that gets no haplotype votes: it is
# left untagged, its confidence column (4) is "." and its PQ column (10) is 0.
# A second read carries the ALT allele of both phased SNPs and must still be
# tagged with a numeric confidence. Needs samtools to build the indexed BAM and
# is skipped without it.
# Usage: tests/test_haplotag_log_no_votes.sh   (after `make`; LONGPHASE overrides the binary)
set -uo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${LONGPHASE:-$root/longphase-to}"
if ! command -v samtools >/dev/null 2>&1; then
    echo "SKIP: haplotag log without votes (samtools not found)"
    exit 0
fi
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cd "$work" || exit 1

ref="$(printf 'ACGTTGCAGTCA%.0s' $(seq 25))"
printf '>chr1\n%s\n' "$ref" > ref.fa
printf 'chr1\t300\t6\t300\t301\n' > ref.fa.fai

# phased SNPs at 100 and 200 (1-based), 0|1 in one phase set
read_seq="${ref:89:120}"
{
    printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=300>\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    for pos in 100 200; do
        base="${ref:pos-1:1}"
        if [ "$base" = G ]; then alt=T; else alt=G; fi
        printf 'chr1\t%s\t.\t%s\t%s\t30\tPASS\t.\tGT:PS\t0|1:100\n' "$pos" "$base" "$alt"
        offset=$((pos - 90))
        read_seq="${read_seq:0:offset}${alt}${read_seq:offset+1}"
    done
} > phased.vcf

# no_votes: 10-59, ends before the first SNP; alt_both: 90-209, ALT at both SNPs
{
    printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n'
    printf 'no_votes\t0\tchr1\t10\t60\t50M\t*\t0\t0\t%s\t*\n' "${ref:9:50}"
    printf 'alt_both\t0\tchr1\t90\t60\t120M\t*\t0\t0\t%s\t*\n' "$read_seq"
} | samtools view -b -o in.bam - && samtools index in.bam

timeout 60 "$bin" haplotag -s phased.vcf -r ref.fa -b in.bam --log -o out > haplotag.log 2>&1
status=$?
if [ "$status" -ne 0 ]; then
    echo "FAIL: haplotag exited with status $status"
    tail -n 3 haplotag.log | sed 's/^/    /'
    exit 1
fi

fail=0
# check <read> <column> <expected>
check() {
    local got
    got="$(awk -F'\t' -v r="$1" -v c="$2" '$1 == r {print $c}' out.out)"
    if [ "$got" != "$3" ]; then
        echo "FAIL: $1 column $2: expected '$3', got '$got'"
        fail=1
    fi
}
check no_votes 4 "."
check no_votes 5 "."
check no_votes 10 "0"
check alt_both 4 "1"
check alt_both 5 "2"
if [ "$fail" -eq 0 ]; then
    echo "PASS: haplotag log without votes"
fi
exit "$fail"
