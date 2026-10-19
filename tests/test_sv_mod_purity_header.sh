#!/usr/bin/env bash
# Check that phase writes the same ##tumor_purity into <prefix>_SV.vcf and
# <prefix>_mod.vcf as into <prefix>.vcf. Uses a two-SNP VCF, header-only SV and
# MOD VCFs and a BAM without reads, with --purity 0.5. Needs samtools to build
# the indexed BAM and is skipped without it.
# Usage: tests/test_sv_mod_purity_header.sh   (after `make`; LONGPHASE overrides the binary)
set -uo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${LONGPHASE:-$root/longphase-to}"
if ! command -v samtools >/dev/null 2>&1; then
    echo "SKIP: SV/MOD purity header (samtools not found)"
    exit 0
fi
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cd "$work" || exit 1

ref="$(printf 'ACGTTGCAGTCA%.0s' $(seq 25))"
printf '>chr1\n%s\n' "$ref" > ref.fa
printf 'chr1\t300\t6\t300\t301\n' > ref.fa.fai
header='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=300>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
{
    printf "$header"
    printf '##source=ClairS-TO\n'
    printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
    printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    printf '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n'
    printf '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    for pos in 100 200; do
        base="${ref:pos-1:1}"
        if [ "$base" = G ]; then alt=T; else alt=G; fi
        printf 'chr1\t%s\t.\t%s\t%s\t30\tPASS\tH\tGT:GQ:DP:AF:AD\t0/1:30:20:0.5:10,10\n' "$pos" "$base" "$alt"
    done
} > snv.vcf
{ printf "$header"; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'; } > sv.vcf
cp sv.vcf mod.vcf
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n' | samtools view -b -o in.bam - && samtools index in.bam

timeout 60 "$bin" phase -s snv.vcf -r ref.fa -b in.bam -c clairs_to_ss --purity 0.5 \
    --sv-file sv.vcf --mod-file mod.vcf -o out > phase.log 2>&1
status=$?
fail=0
if [ "$status" -ne 0 ]; then
    echo "FAIL: phase exited with status $status"
    tail -n 3 phase.log | sed 's/^/    /'
    exit 1
fi
for f in out.vcf out_SV.vcf out_mod.vcf; do
    value="$(grep -m1 '^##tumor_purity=' "$f" 2>/dev/null)"
    if [ "$value" != "##tumor_purity=0.5" ]; then
        echo "FAIL: $f: expected ##tumor_purity=0.5, got '${value}'"
        fail=1
    fi
done
if [ "$fail" -eq 0 ]; then
    echo "PASS: SV/MOD purity header"
fi
exit "$fail"
