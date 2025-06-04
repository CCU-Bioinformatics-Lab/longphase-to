## Input Preparation
- [Generate reference index](./docs/input_preparation.md#generate-reference-index)
- [Generate alignment and index files](./docs/input_preparation.md#generate-alignment-and-index-files)
- [Generate single nucleotide polymorphism (SNP) file](#generate-single-nucleotide-polymorphism-snp-file)
    - [ClairS-TO Caller](#clairs-to-caller)
    - [DeepSomatic Caller](#deepsomatic-caller)
- [Downloads panels of normals (PoNs) file](#downloads-panels-of-normals-pons-file)
- [Generate Structural variation (SV) file](./docs/input_preparation.md#generate-structural-variation-sv-file)
- [Carry methylation tags to BAMs](./docs/input_preparation.md#carry-methylation-tags-to-bams)
### Generate reference index
Index the reference genome with [samtools](https://github.com/samtools/samtools).
```
samtools faidx reference.fasta
```
### Generate alignment and index files
Produce read-to-reference alignment via [minimap2](https://github.com/lh3/minimap2) and sort/index the bam by [samtools](https://github.com/samtools/samtools).
```
# generate alignment flie with minimap2 according to the sequencing platform e.g. map-pb/map-ont/map-hifi
# Note that the MD-tag is required by sniffles (–MD).
minimap2 --MD -ax map-ont -t 10 reference.fasta reads.fastq -o alignment.sam

# sort alignment file
samtools sort -@ 10 alignment.sam -o alignment.bam

# index alignment file
samtools index -@ 10 alignment.bam
```
### Generate single nucleotide polymorphism (SNP) file
#### ClairS-TO Caller
```
INPUT_BAM_DIR="/path/to/bam"
INPUT_REF_DIR="/path/to/reference"
OUTPUT_DIR="/path/to/output"
BAM="alignment.bam"
REF="reference.fasta"
THREADS=64
MODEL="ont_r10_dorado_sup_5khz_ssrs"
# ont_r10_dorado_sup_5khz_ssrs is preset for ONT R10 Dorado 5 "Sup" basecaller
# ssrs is a model trained initially with synthetic samples and then real samples augmented
# if you do not want to use real data, use the ss model.
# for PacBio-HiFi reads: hifi_revio_ssrs

sudo docker run \
-v ${INPUT_BAM_DIR}:${INPUT_BAM_DIR} \
-v ${INPUT_REF_DIR}:${INPUT_REF_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
-u $(id -u):$(id -g) \
hkubal/clairs-to:v0.3.0 \
/opt/bin/run_clairs_to \
--tumor_bam_fn ${INPUT_BAM_DIR}/${BAM} \
--ref_fn ${INPUT_REF_DIR}/${REF} \
--threads ${THREADS} \
--platform ${MODEL} \
--output_dir ${OUTPUT_DIR}
```

#### DeepSomatic Caller
```
INPUT_BAM_DIR="/path/to/bam"
INPUT_REF_DIR="/path/to/reference"
OUTPUT_DIR="/path/to/output"
BAM="alignment.bam"
REF="reference.fasta"
THREADS=64
MODEL="ONT_TUMOR_ONLY"
# for PacBio-HiFi reads: PACBIO_TUMOR_ONLY

sudo docker run \
-v ${INPUT_BAM_DIR}:${INPUT_BAM_DIR} \
-v ${INPUT_REF_DIR}:${INPUT_REF_DIR} \
-u $(id -u):$(id -g) \
google/deepsomatic:1.8.0 \
run_deepsomatic \
--model_type ${MODEL} \
--ref ${INPUT_REF_DIR}/${REF} \
--reads_tumor  ${INPUT_BAM_DIR}/${BAM} \
--output_vcf  ${OUTPUT_DIR}/output.vcf.gz \
--sample_name_tumor "tumor" \
--num_shards ${THREADS} \
--logging_dir ${OUTPUT_DIR}/logs \
--intermediate_results_dir ${OUTPUT_DIR}/intermediate_results_dir \
--use_default_pon_filtering=true
```

If you're using GPUs.
```
INPUT_BAM_DIR="/path/to/bam"
INPUT_REF_DIR="/path/to/reference"
OUTPUT_DIR="/path/to/output"
BAM="alignment.bam"
REF="reference.fasta"
THREADS=64
MODEL="ONT_TUMOR_ONLY"
# for PacBio-HiFi reads: PACBIO_TUMOR_ONLY

sudo docker run --gpus all \
-v ${INPUT_BAM_DIR}:${INPUT_BAM_DIR} \
-v ${INPUT_REF_DIR}:${INPUT_REF_DIR} \
-u $(id -u):$(id -g) \
google/deepsomatic:1.8.0-gpu \
run_deepsomatic \
--model_type ${MODEL} \
--ref ${INPUT_REF_DIR}/${REF} \
--reads_tumor  ${INPUT_BAM_DIR}/${BAM} \
--output_vcf  ${OUTPUT_DIR}/output.vcf.gz \
--sample_name_tumor "tumor" \
--num_shards ${THREADS} \
--logging_dir ${OUTPUT_DIR}/logs \
--intermediate_results_dir ${OUTPUT_DIR}/intermediate_results_dir \
--use_default_pon_filtering=true
```

### Generate panels of normals (PoNs) file
source: [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO?tab=readme-ov-file#tagging-non-somatic-variant-using-panels-of-normals-pons)
```
PATH="PoN"
mkdir $PATH
wget -P $PATH http://www.bio8.cs.hku.hk/clairs-to/databases/gnomad.r2.1.af-ge-0.001.sites.vcf.gz
wget -P $PATH http://www.bio8.cs.hku.hk/clairs-to/databases/dbsnp.b138.non-somatic.sites.vcf.gz
wget -P $PATH http://www.bio8.cs.hku.hk/clairs-to/databases/1000g-pon.sites.vcf.gz
wget -P $PATH http://www.bio8.cs.hku.hk/clairs-to/databases/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz

# longphase-to phase \
# --pon-file PoN/1000g-pon.sites.vcf.gz,\
#            PoN/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz,\
# --strict-pon-file PoN/dbsnp.b138.non-somatic.sites.vcf.gz,\
#                   PoN/gnomad.r2.1.af-ge-0.001.sites.vcf.gz
```

### Generate Structural variation (SV) file
e.g. [sniffles](https://github.com/fritzsedlazeck/Sniffles) or [CuteSV](https://github.com/tjiangHIT/cuteSV).
```
# In sniffles1 please specofic --num_reads_report -1. For sniffles2 please specify --output-rnames instead.
sniffles -t 10 --num_reads_report -1 -m alignment.bam -v SV.vcf # for sniffles1
sniffles --threads 10 --output-rnames --input alignment.bam --vcf SV.vcf # for sniffles2

# cuteSV command for PacBio CLR data:
cuteSV alignment.bam reference.fasta SV.vcf work_dir --report_readid --genotype

# additional platform-specific parameters suggested by cuteSV
# PacBio CLR data: 
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5
# PacBio CCS(HIFI) data: 
--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
# ONT data: 
--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3
```

#### Carry methylation tags to BAMs
[The `-T` parameter](https://github.com/nanoporetech/dorado/issues/145) in samtools fastq extracts tags from the BAM file and stores them in the header of the FASTQ file. Please ensure that the BAM file includes both `MM` and `ML` tags and carried on in the following way.
```
samtools fastq -T '*' methylcall.raw.bam > methylcall.raw.fastq
```

Then, specify the `-y` option in minimap2 which appends tags stored in the FASTQ header into the BAM file.
```
minimap2 -ax map-ont -y reference.fasta methylcall.raw.fastq 
```
