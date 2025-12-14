#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------------------------
# Cohesive pipeline for targeted extraction & variant calling
#   - Validates inputs and tools
#   - Harmonizes contig naming ("chr"/no-"chr") per BAM
#   - Extracts reads overlapping a region
#   - Optional local re-alignment to the reference
#   - bcftools mpileup -> call (-c by default) per sample
#   - BGZF-compressed, indexed per-sample VCFs
#   - Multi-sample merge and genotype matrix (GT; optional DP)
#   - Optional single-site (pos) per-sample summary (GT, DP, depth fallback)
#
# Requirements: bash (>=4.3 for wait -n), samtools, bcftools, bwa (if --realign),
#               bgzip/tabix (or bcftools index), coreutils
# -------------------------------------------------------------------

# Defaults
THREADS=4
CALLER_MODE="c"          # bcftools call mode: "c" (legacy consensus) or "mv"
REALIGN="no"             # "yes" to realign the regional BAMs with BWA-MEM
INPUT_LIST=""            # file with absolute/relative paths to BAMs (one per line)
INPUT_DIR=""             # alternatively, directory to scan for *.bam
REF_FA=""                # reference FASTA (requires .fai; if realigning, requires BWA index)
REGION=""                # e.g. "chr20:5904018-5904038"
SINGLE_SITE=""           # optional e.g. "chr20:5904028"
OUTDIR="out_target"
KEEP_INTERMEDIATE="no"   # "yes" keeps regional/realigned BAMs
GT_ONLY="yes"            # "no" -> also output a DP matrix
LOGDIR="logs"

usage() {
  cat <<EOF
Usage:
  $(basename "$0") -r REF.fa -g REGION [-s SINGLE_SITE] (-i bam_list.txt | -d bam_dir) [options]

Required:
  -r  Reference FASTA (indexed with samtools faidx; if --realign, also bwa-indexed)
  -g  Target region (e.g. chr20:5904018-5904038)
  -i  Text file with BAM paths (one per line)   OR
  -d  Directory to scan for *.bam

Optional:
  -s  Single genomic position for summary (e.g. chr20:5904028)
  -o  Output directory (default: ${OUTDIR})
  -t  Threads (default: ${THREADS})
  -m  Caller mode: c (default) or mv  [bcftools call -c  OR  -m -v]
  --realign            Realign extracted reads with BWA-MEM (default: no)
  --keep-intermediate  Keep regional/realigned BAMs (default: no)
  --gt-only=no         Also export a DP matrix (default: GT-only)
  -h  Show this help

Examples:
  $(basename "$0") -r ref.fa -g chr20:5904018-5904038 -s chr20:5904028 -i bam_list.txt -t 8 --realign
  $(basename "$0") -r ref.fa -g 20:100000-101000 -d /data/bams -m mv --gt-only=no
EOF
  exit 1
}

# Parse args
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r) REF_FA="$2"; shift 2;;
    -g) REGION="$2"; shift 2;;
    -s) SINGLE_SITE="$2"; shift 2;;
    -i) INPUT_LIST="$2"; shift 2;;
    -d) INPUT_DIR="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -m) CALLER_MODE="$2"; shift 2;;
    --realign) REALIGN="yes"; shift 1;;
    --keep-intermediate) KEEP_INTERMEDIATE="yes"; shift 1;;
    --gt-only=no) GT_ONLY="no"; shift 1;;
    -h|--help) usage;;
    *) ARGS+=("$1"); shift 1;;
  esac
done

# Validate required inputs
[[ -z "${REF_FA}" || -z "${REGION}" ]] && usage
[[ -z "${INPUT_LIST}" && -z "${INPUT_DIR}" ]] && usage
[[ -n "${INPUT_LIST}" && -n "${INPUT_DIR}" ]] && { echo "Provide either -i or -d, not both." >&2; exit 1; }

# Tool checks
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing required tool: $1" >&2; exit 1; }; }
need samtools
need bcftools
if [[ "${REALIGN}" == "yes" ]]; then need bwa; fi

# Reference checks
[[ -f "${REF_FA}" ]] || { echo "Reference FASTA not found: ${REF_FA}" >&2; exit 1; }
[[ -f "${REF_FA}.fai" ]] || { echo "Missing FASTA index: ${REF_FA}.fai (run: samtools faidx ${REF_FA})" >&2; exit 1; }
if [[ "${REALIGN}" == "yes" ]]; then
  for ext in amb ann bwt pac sa; do
    [[ -f "${REF_FA}.${ext}" ]] || { echo "Missing BWA index ${REF_FA}.${ext} (run: bwa index ${REF_FA})" >&2; exit 1; }
  done
fi

# Prepare I/O
mkdir -p "${OUTDIR}" "${LOGDIR}"
SAMPLES_DIR="${OUTDIR}/samples"
mkdir -p "${SAMPLES_DIR}"

# Collect BAM list
BAMS_TXT="${OUTDIR}/bam_list.resolved.txt"
> "${BAMS_TXT}"
if [[ -n "${INPUT_LIST}" ]]; then
  while IFS='' read -r line; do
    [[ -z "${line}" ]] && continue
    [[ "${line}" =~ ^# ]] && continue
    echo "${line}"
  done < "${INPUT_LIST}" > "${BAMS_TXT}"
else
  find "${INPUT_DIR}" -type f -name "*.bam" -print0 | xargs -0 -I{} echo "{}" >> "${BAMS_TXT}"
fi

# Sanity check BAMs
echo "[INFO] Verifying BAM paths and quick integrity..." >&2
VALID_BAMS="${OUTDIR}/bam_list.valid.txt"
> "${VALID_BAMS}"
while IFS='' read -r BAM; do
  if [[ ! -f "${BAM}" ]]; then
    echo "[WARN] Missing: ${BAM}" >&2
    continue
  fi
  if ! samtools quickcheck -v "${BAM}" 2>>"${LOGDIR}/quickcheck.err"; then
    echo "[WARN] samtools quickcheck failed: ${BAM}" >&2
    continue
  fi
  echo "${BAM}" >> "${VALID_BAMS}"
done < "${BAMS_TXT}"

# Helpers
bam_has_chr_prefix() {  # returns 0 if header SN contain 'chr' prefix
  samtools view -H "$1" | grep -q "^@SQ.*\sSN:chr"
}

adapt_region_for_bam() { # $1=BAM, $2=region -> echo adjusted region
  local bam="$1"; local reg="$2"
  local chrom="${reg%%:*}"
  local rest="${reg#*:}"; [[ "${rest}" == "${reg}" ]] && rest=""
  local has_chr_header=1
  if bam_has_chr_prefix "${bam}"; then has_chr_header=0; fi

  # Determine if region has chr
  local chrom_has_chr=1
  [[ "${chrom}" == chr* ]] && chrom_has_chr=0

  if [[ ${has_chr_header} -eq 0 && ${chrom_has_chr} -ne 0 ]]; then
    chrom="chr${chrom}"
  elif [[ ${has_chr_header} -ne 0 && ${chrom_has_chr} -eq 0 ]]; then
    chrom="${chrom#chr}"
  fi

  if [[ -n "${rest}" && "${rest}" != "${reg}" ]]; then
    echo "${chrom}:${rest}"
  else
    echo "${chrom}"
  fi
}

ensure_bam_index() { # $1=BAM
  local b="$1"
  if [[ ! -f "${b}.bai" && ! -f "${b}.csi" ]]; then
    samtools index -@ "${THREADS}" "${b}"
  fi
}

sample_name_from_bam() { # $1=BAM -> echo sample name
  local base
  base="$(basename "$1")"
  echo "${base%.bam}"
}

# Variant calling command builder
build_call_cmd() {
  local mode="$1"
  if [[ "${mode}" == "mv" ]]; then
    echo "bcftools call -m -v"
  else
    echo "bcftools call -c"
  fi
}

# Per-sample processing
process_bam() {
  local BAM="$1"
  local SAMPLE; SAMPLE="$(sample_name_from_bam "${BAM}")"
  local SDIR="${SAMPLES_DIR}/${SAMPLE}"
  mkdir -p "${SDIR}"

  local LOG="${LOGDIR}/${SAMPLE}.log"
  echo "[INFO] Processing ${SAMPLE}" | tee "${LOG}"

  ensure_bam_index "${BAM}"

  local ADJ_REGION; ADJ_REGION="$(adapt_region_for_bam "${BAM}" "${REGION}")"
  echo "[INFO] Adjusted region for ${SAMPLE}: ${ADJ_REGION}" | tee -a "${LOG}"

  # Extract regional BAM
  local REG_BAM="${SDIR}/${SAMPLE}.region.bam"
  samtools view -@ "${THREADS}" -b -o "${REG_BAM}" "${BAM}" "${ADJ_REGION}" >> "${LOG}" 2>&1 || true
  # If empty, we still proceed to create an index (may be empty); note of zero reads.
  if [[ ! -s "${REG_BAM}" ]]; then
    echo "[WARN] No overlapping reads for ${SAMPLE} in ${ADJ_REGION}" | tee -a "${LOG}"
  fi
  samtools index -@ "${THREADS}" "${REG_BAM}" >> "${LOG}" 2>&1 || true

  # Optional local realignment
  local ALN_BAM="${REG_BAM}"
  if [[ "${REALIGN}" == "yes" && -s "${REG_BAM}" ]]; then
    ALN_BAM="${SDIR}/${SAMPLE}.realigned.bam"
    samtools fastq -@ "${THREADS}" "${REG_BAM}" \
    | bwa mem -t "${THREADS}" "${REF_FA}" - \
    | samtools sort -@ "${THREADS}" -o "${ALN_BAM}" >> "${LOG}" 2>&1
    samtools index -@ "${THREADS}" "${ALN_BAM}" >> "${LOG}" 2>&1
  fi

  # Variant calling (skip if no data)
  local VCF_GZ="${SDIR}/${SAMPLE}.vcf.gz"
  if [[ -s "${ALN_BAM}" ]]; then
    local ADJ_REGION_FOR_ALN; ADJ_REGION_FOR_ALN="$(adapt_region_for_bam "${ALN_BAM}" "${REGION}")"
    local CALL_CMD; CALL_CMD="$(build_call_cmd "${CALLER_MODE}")"
    # mpileup -> call -> bgzip
    bcftools mpileup -f "${REF_FA}" -a AD,DP -r "${ADJ_REGION_FOR_ALN}" -Ou "${ALN_BAM}" \
      | eval "${CALL_CMD}" -Ou \
      | bcftools view -Oz -o "${VCF_GZ}" >> "${LOG}" 2>&1 || true

    if [[ -s "${VCF_GZ}" ]]; then
      bcftools index -t -f "${VCF_GZ}" >> "${LOG}" 2>&1 || true
    else
      echo "[WARN] Empty VCF for ${SAMPLE}" | tee -a "${LOG}"
      rm -f "${VCF_GZ}" "${VCF_GZ}.tbi" || true
    fi
  else
    echo "[WARN] Skipping call; empty alignment for ${SAMPLE}" | tee -a "${LOG}"
  fi

  # Optional cleanup
  if [[ "${KEEP_INTERMEDIATE}" != "yes" ]]; then
    rm -f "${REG_BAM}" "${REG_BAM}.bai" || true
    if [[ "${REALIGN}" == "yes" ]]; then
      rm -f "${ALN_BAM}" "${ALN_BAM}.bai" || true
    fi
  fi
}

# Parallel per-sample loop
run_per_sample() {
  local njobs=0
  while IFS='' read -r BAM; do
    [[ -z "${BAM}" ]] && continue
    process_bam "${BAM}" &
    ((njobs++))
    if (( njobs >= THREADS )); then
      wait -n
      ((njobs--))
    fi
  done < "${VALID_BAMS}"
  wait
}

# Merge per-sample VCFs (only those produced)
merge_vcfs() {
  local MERGED="${OUTDIR}/merged.vcf.gz"
  local LIST="${OUTDIR}/vcf_list.txt"
  find "${SAMPLES_DIR}" -type f -name "*.vcf.gz" | sort > "${LIST}" || true

  if [[ ! -s "${LIST}" ]]; then
    echo "[WARN] No per-sample VCFs to merge." >&2
    return 0
  fi

  # If only one VCF, copy; else merge
  local COUNT
  COUNT=$(wc -l < "${LIST}")
  if (( COUNT == 1 )); then
    cp "$(head -n1 "${LIST}")" "${MERGED}"
    bcftools index -t -f "${MERGED}" >/dev/null 2>&1 || true
  else
    bcftools merge -m all -Oz -o "${MERGED}" -l "${LIST}"
    bcftools index -t -f "${MERGED}" >/dev/null 2>&1 || true
  fi
  echo "${MERGED}"
}

# Export matrices
export_matrices() {
  local MERGED="$1"
  if [[ -z "${MERGED}" || ! -f "${MERGED}" ]]; then
    echo "[WARN] Skipping matrix export; merged VCF missing." >&2
    return 0
  fi
  local GT_TSV="${OUTDIR}/matrix.genotypes.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "${MERGED}" > "${GT_TSV}"

  if [[ "${GT_ONLY}" == "no" ]]; then
    local DP_TSV="${OUTDIR}/matrix.depth.tsv"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' "${MERGED}" > "${DP_TSV}"
  fi
}

# Single-site per-sample summary: Sample,GT,DP,(DepthFallback)
single_site_summary() {
  local SITE="$1"
  [[ -z "${SITE}" ]] && return 0
  local CSV="${OUTDIR}/single_site_summary.csv"
  echo "Sample,GT,DP,DepthFallback" > "${CSV}"

  while IFS='' read -r BAM; do
    [[ -z "${BAM}" ]] && continue
    local SAMPLE; SAMPLE="$(sample_name_from_bam "${BAM}")"
    local SDIR="${SAMPLES_DIR}/${SAMPLE}"
    local VCF_GZ="${SDIR}/${SAMPLE}.vcf.gz"
    local GT="NA"; local DP="NA"; local DEPTHFB="NA"

    if [[ -f "${VCF_GZ}" ]]; then
      # Limit to this sample; if empty, keep NA
      local Q
      Q=$(bcftools query -s "${SAMPLE}" -r "${SITE}" -f '[%GT\t%DP]\n' "${VCF_GZ}" || true)
      if [[ -n "${Q}" ]]; then
        GT="$(echo "${Q}" | awk '{print $1}')"
        DP="$(echo "${Q}" | awk '{print $2}')"
      fi
    fi

    # If GT still NA, provide a depth fallback from BAM (regional or realigned may be removed)
    # Use original BAM for robustness.
    local ADJ_SITE; ADJ_SITE="$(adapt_region_for_bam "${BAM}" "${SITE}")"
    local DEP
    DEP=$(samtools depth -r "${ADJ_SITE}" "${BAM}" | awk '{sum=$3} END{if (sum=="") print 0; else print sum}')
    DEPTHFB="${DEP}"

    echo "${SAMPLE},${GT},${DP},${DEPTHFB}" >> "${CSV}"
  done < "${VALID_BAMS}"
}

# -------------------------
# Run the pipeline
# -------------------------
echo "[INFO] Starting per-sample processing..." >&2
run_per_sample

echo "[INFO] Merging per-sample VCFs..." >&2
MERGED_VCF="$(merge_vcfs || true)"
[[ -n "${MERGED_VCF:-}" ]] && echo "[INFO] Merged VCF: ${MERGED_VCF}" >&2

echo "[INFO] Exporting matrices..." >&2
export_matrices "${MERGED_VCF:-}"

if [[ -n "${SINGLE_SITE}" ]]; then
  echo "[INFO] Building single-site summary for ${SINGLE_SITE}..." >&2
  single_site_summary "${SINGLE_SITE}"
fi

echo "[INFO] Done. Outputs in: ${OUTDIR}" >&2
