#!/usr/bin/env bash
# Check that phase and haplotag stop with an error message, not an abort, when a
# contig they need is missing from the reference. phase needs the contigs that
# carry SNPs; haplotag needs every contig in the phased VCF header. The
# reference here has chr1 only. The haplotag check needs samtools to build an
# indexed BAM and is skipped without it.
# Usage: tests/test_missing_contig.sh   (after `make`; LONGPHASE overrides the binary)
set -uo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${LONGPHASE:-$root/longphase-to}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cd "$work" || exit 1
ulimit -c 0

printf '>chr1\nACGTACGTACGTACGTACGT\n' > ref.fa
printf 'chr1\t20\t6\t20\t21\n' > ref.fa.fai

vcf_header() {
    printf '##fileformat=VCFv4.2\n##source=ClairS-TO\n'
    printf '##contig=<ID=chr1,length=20>\n##contig=<ID=chr2,length=20>\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
    printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    printf '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n'
    printf '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
    printf '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
}
fail=0

# expect_error <name> <command...>: exit status 1 and a message naming chr2
expect_error() {
    local name="$1" err status
    shift
    err="$(timeout 60 "$@" 2>&1 >/dev/null)"
    status=$?
    if [ "$status" -ne 1 ] || ! grep -q "chr2 not found in reference" <<<"$err"; then
        echo "FAIL: $name: expected exit status 1 and 'chr2 not found in reference', got status $status"
        tail -n 3 <<<"$err" | sed 's/^/    /'
        fail=1
    else
        echo "PASS: $name"
    fi
}

# phase: SNP on chr2, which the reference lacks
{ vcf_header; printf 'chr2\t5\t.\tA\tG\t30\tPASS\tH\tGT:GQ:DP:AF:AD\t0/1:30:20:0.5:10,10\n'; } > snv.vcf
: > in.bam
expect_error "phase, SNP on a contig missing from the reference" \
    "$bin" phase -s snv.vcf -r ref.fa -b in.bam -c clairs_to_ss --purity 0.5 -o phase_out

# haplotag: chr2 is only in the VCF header
if command -v samtools >/dev/null 2>&1; then
    { vcf_header; printf 'chr1\t5\t.\tA\tG\t30\tPASS\tH\tGT:GQ:DP:AF:AD:PS\t0|1:30:20:0.5:10,10:5\n'; } > phased.vcf
    printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n' | samtools view -b -o tag.bam - && samtools index tag.bam
    expect_error "haplotag, VCF header contig missing from the reference" \
        "$bin" haplotag -s phased.vcf -r ref.fa -b tag.bam -o tag_out
else
    echo "SKIP: haplotag check (samtools not found)"
fi

if [ "$fail" -eq 0 ]; then
    echo "PASS: missing contig"
fi
exit "$fail"
