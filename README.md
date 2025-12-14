Targeted WES Extractor & Splice Junction Pipeline

This repository contains two complementary tools for extracting targeted genomic regions from WES data and quantifying splice junction usage from aligned RNA/DNA sequencing BAM files.

Contents

WesExtractor.sh — a Bash pipeline for targeted extraction, re-alignment, and variant calling within user-defined genomic windows.

SpliceJunctionFromBam.py — a Python utility for detecting splice junctions directly from BAM files using CIGAR parsing, outputting regtools-compatible BED12 junction tables.

1. WesExtractor.sh

WesExtractor.sh performs automated extraction of a genomic locus from multiple WES BAM files, standardizing coordinates, optionally re-aligning extracted reads, and generating per-sample VCF files suitable for downstream statistical or clinical workflows.

Features

Validates BAM paths and indexes if missing

Extracts reads overlapping a user-defined region (e.g., chr20:5904018-5904038)

Optional BWA-MEM realignment of regional reads

Per-sample variant calling via bcftools mpileup → call

Compressed, indexed VCF output for rapid querying

Summary genotype/coverage extraction for key loci

Example
bash WesExtractor.sh -i bam_list.txt \
                     -r reference.fa \
                     -g chr20:5904018-5904038 \
                     -t 8 \
                     --realign

Output

sample.region.bam (optional)

sample.vcf.gz and index

summary.csv with GT and DP per position

2. SpliceJunctionFromBam.py

SpliceJunctionFromBam.py extracts splice junctions by scanning CIGAR strings for “N” operations and enforcing biologically relevant anchor and intron-length thresholds. The script outputs BED12 formatted junctions compatible with regtools and downstream quantitative splicing tools.
The logic is identical to that used in the unified MAPT junction module.

Features

Uses pysam to iterate over reads in a target genomic region

Detects splice junctions using CIGAR parsing

Applies configurable constraints:

Minimum anchor length

Intron-length bounds

Produces *.MAPT.juncs.bed (or general BED12) files

Example
python3 SpliceJunctionFromBam.py \
    --bam sample.bam \
    --region 17:45890000-46050000 \
    --out sample.MAPT.juncs.bed

Output Format

A BED12 file with:

Two exon blocks representing the splice donor and acceptor

Junction read count in the score field

Coordinates compatible with regtools and junction-quantification workflows

Requirements
For WesExtractor.sh

samtools, bcftools, bgzip/tabix, bwa (if re-alignment enabled)

POSIX-compatible shell environment (Linux, WSL, or macOS)

For SpliceJunctionFromBam.py

Python ≥ 3.8

pysam
