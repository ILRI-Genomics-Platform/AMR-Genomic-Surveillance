#!/usr/bin/env bash
# run_mpox_steps3-10.sh
# Automate Steps 3 -> 10 from mpox-genome-analysis-illumina.md
# Usage:
#   ./run_mpox_steps3-10.sh              # autodetect samples from DATA_DIR
#   ./run_mpox_steps3-10.sh SRR12345 SRR67890  # run only specified sample IDs
#
# IMPORTANT:
# - Load modules / activate conda env as in Step 1/2 of the tutorial before running.
# - Edit the configuration block below to suit your filesystem and resource limits.

set -euo pipefail
IFS=$'\n\t'

### -----------------------
### Configuration (edit)
### -----------------------
# Base directories (relative to the repo root or absolute)
DATA_DIR="./data/mpox/illumina"                  # where raw fastq are stored
RESULTS_DIR="./results/mpox"                     # base results folder
PRIMER_SCHEME_DIR="./primer_scheme"              # where reference.fasta and primer.bed live (source)
PRIMER_SCHEME_OUT="${RESULTS_DIR}/primerschemes" # copied scheme location used by pipeline
BWA_INDEX_DIR="${RESULTS_DIR}/bwa/index"
BWA_ALIGN_DIR="${RESULTS_DIR}/bwa/alignment"
HOSTILE_INDEX="./databases/hostile/human-t2t-hla" # hostile index path
THREADS=4
SAMTOOLS_THREADS=2
NORMALISE=200
PRIMER_MATCH_THRESHOLD=35
MIN_MAPQ=20
FREEBAYES_MIN_COVERAGE=10
FREEBAYES_FREQ_LOW=0.2
FREEBAYES_MIN_COUNT=1

# Squirrel / step10 settings
RUN_SQUIRREL=false
SQUIRREL_THREADS=1

# Safety: commands used (script will check they exist)
REQUIRED_CMDS=(hostile fastp bwa samtools freebayes bgzip tabix bcftools squirrel)

### -----------------------
### Helper / preparation
### -----------------------
echolog(){ echo "[$(date -u +'%Y-%m-%dT%H:%M:%SZ')] $*"; }

# check required tools (warn only for optional ones)
for c in "${REQUIRED_CMDS[@]}"; do
    if ! command -v "$c" >/dev/null 2>&1; then
        # allow missing squirrel (optional)
        if [ "$c" = "squirrel" ] && [ "${RUN_SQUIRREL}" = "false" ]; then
            continue
        fi
        echo "ERROR: Required command not found in PATH: $c" >&2
        echo "Please load modules or activate environment containing: ${REQUIRED_CMDS[*]}" >&2
        exit 1
    fi
done

# create required directory tree that matches tutorial
mkdir -p \
  "${RESULTS_DIR}/hostile" \
  "${RESULTS_DIR}/fastp" \
  "${PRIMER_SCHEME_OUT}" \
  "${BWA_INDEX_DIR}" \
  "${BWA_ALIGN_DIR}" \
  "${RESULTS_DIR}/primertrimmed" \
  "${RESULTS_DIR}/freebayes" \
  "${RESULTS_DIR}/data/ncbi" \
  "${RESULTS_DIR}/data/pathoplexus" \
  "${RESULTS_DIR}/data/all-consensus" \
  "${RESULTS_DIR}/squirrel"

### -----------------------
### Sample list detection
### -----------------------
SAMPLES=()
if [ $# -gt 0 ]; then
  # samples supplied on command line
  SAMPLES=( "$@" )
else
  # autodetect by looking for *_1.fastq.gz
  for fq1 in "${DATA_DIR}"/*_1.fastq.gz "${DATA_DIR}"/*_R1_*.fastq.gz 2>/dev/null; do
    [ -e "$fq1" ] || continue
    # handle names like SRRxxx_1.fastq.gz  -> sample prefix before _1.fastq.gz
    base=$(basename "$fq1")
    if [[ "$base" =~ (.*)_1\.fastq\.gz$ ]]; then
      SAMPLES+=( "${BASH_REMATCH[1]}" )
    elif [[ "$base" =~ (.*)_R1_.*\.fastq\.gz$ ]]; then
      # e.g. sample_R1_001.fastq.gz -> sample
      SAMPLES+=( "${BASH_REMATCH[1]}" )
    fi
  done
fi

if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "No samples detected in ${DATA_DIR}. Provide sample IDs as arguments or place fastq files as <SAMPLE>_1.fastq.gz and <SAMPLE>_2.fastq.gz"
  exit 1
fi

echolog "Samples to process: ${SAMPLES[*]}"

### -----------------------
### Step 5 prep: copy primer scheme & index reference if missing
### We copy reference.fasta and primer.bed into PRIMER_SCHEME_OUT to mimic tutorial
### -----------------------
if [ ! -f "${PRIMER_SCHEME_OUT}/reference.fasta" ] || [ ! -f "${PRIMER_SCHEME_OUT}/primer.bed" ]; then
  echolog "Copying primer scheme files to ${PRIMER_SCHEME_OUT}"
  cp -v "${PRIMER_SCHEME_DIR}/reference.fasta" "${PRIMER_SCHEME_OUT}/" || { echolog "ERROR: reference.fasta not found in ${PRIMER_SCHEME_DIR}"; exit 1; }
  cp -v "${PRIMER_SCHEME_DIR}/primer.bed" "${PRIMER_SCHEME_OUT}/" || { echolog "ERROR: primer.bed not found in ${PRIMER_SCHEME_DIR}"; exit 1; }
fi

# Create BWA index if missing
if [ ! -f "${BWA_INDEX_DIR}/reference.amb" ]; then
  echolog "Creating BWA index in ${BWA_INDEX_DIR} (this runs once)"
  mkdir -p "${BWA_INDEX_DIR}"
  bwa index -p "${BWA_INDEX_DIR}/reference" "${PRIMER_SCHEME_OUT}/reference.fasta"
fi

BWA_INDEX="${BWA_INDEX_DIR}/reference"

### -----------------------
### Per-sample processing functions (Steps 3 -> 9)
### -----------------------
process_sample() {
  local SAMPLE="$1"
  echolog "=== START SAMPLE: ${SAMPLE} ==="

  # file names (match tutorial)
  local RAW_R1="${DATA_DIR}/${SAMPLE}_1.fastq.gz"
  local RAW_R2="${DATA_DIR}/${SAMPLE}_2.fastq.gz"
  # also accept _R1_001 naming
  if [ ! -f "$RAW_R1" ] && [ -f "${DATA_DIR}/${SAMPLE}_R1_001.fastq.gz" ]; then
      RAW_R1="${DATA_DIR}/${SAMPLE}_R1_001.fastq.gz"
      RAW_R2="${DATA_DIR}/${SAMPLE}_R2_001.fastq.gz"
  fi

  [ -f "$RAW_R1" ] || { echolog "Missing R1: $RAW_R1"; return 1; }
  [ -f "$RAW_R2" ] || { echolog "Missing R2: $RAW_R2"; return 1; }

  # Step 3: Remove human reads -> hostile
  echolog "Step 3: removing human reads with hostile for ${SAMPLE}"
  hostile clean \
    --fastq1 "${RAW_R1}" \
    --fastq2 "${RAW_R2}" \
    --output "${RESULTS_DIR}/hostile/" \
    --threads "${THREADS}" \
    --index "${HOSTILE_INDEX}"

  # hostile naming: tutorial uses results/mpox/hostile/SRR21755837_1.clean_1.fastq.gz etc.
  HOST_CLEAN_R1="${RESULTS_DIR}/hostile/${SAMPLE}_1.clean_1.fastq.gz"
  HOST_CLEAN_R2="${RESULTS_DIR}/hostile/${SAMPLE}_2.clean_2.fastq.gz"
  if [ ! -f "${HOST_CLEAN_R1}" ]; then
    # alternate pattern used in some hostile versions
    HOST_CLEAN_R1="${RESULTS_DIR}/hostile/${SAMPLE}_1.clean.fastq.gz"
    HOST_CLEAN_R2="${RESULTS_DIR}/hostile/${SAMPLE}_2.clean.fastq.gz"
  fi

  [ -f "${HOST_CLEAN_R1}" ] || { echolog "Hostile output R1 missing for ${SAMPLE}: expected ${HOST_CLEAN_R1}"; return 1; }
  [ -f "${HOST_CLEAN_R2}" ] || { echolog "Hostile output R2 missing for ${SAMPLE}: expected ${HOST_CLEAN_R2}"; return 1; }

  # Step 4: Trim adapters with fastp
  echolog "Step 4: trimming adapters with fastp for ${SAMPLE}"
  mkdir -p "${RESULTS_DIR}/fastp"
  FASTP_R1="${RESULTS_DIR}/fastp/${SAMPLE}_R1.trim.fastq.gz"
  FASTP_R2="${RESULTS_DIR}/fastp/${SAMPLE}_R2.trim.fastq.gz"
  FASTP_JSON="${RESULTS_DIR}/fastp/${SAMPLE}.fastp.json"
  FASTP_HTML="${RESULTS_DIR}/fastp/${SAMPLE}.fastp.html"
  FASTP_LOG="${RESULTS_DIR}/fastp/${SAMPLE}.fastp.log"

  fastp \
    --in1 "${HOST_CLEAN_R1}" \
    --in2 "${HOST_CLEAN_R2}" \
    --out1 "${FASTP_R1}" \
    --out2 "${FASTP_R2}" \
    --detect_adapter_for_pe \
    --thread "${THREADS}" \
    --json "${FASTP_JSON}" \
    --html "${FASTP_HTML}" \
    --cut_mean_quality 20 \
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 40 \
    --length_required 20 \
    2> "${FASTP_LOG}"

  # Step 5 already handled (we copied scheme & reference earlier)

  # Step 6: Alignment
  echolog "Step 6: aligning reads for ${SAMPLE}"
  mkdir -p "${BWA_ALIGN_DIR}"
  # Align to reference -> produce BAM
  BWA_BAM="${BWA_ALIGN_DIR}/${SAMPLE}.bam"
  SAMTOOLS_SORTED="${BWA_ALIGN_DIR}/${SAMPLE}.sorted.bam"

  bwa mem -t "${THREADS}" "${BWA_INDEX}" "${FASTP_R1}" "${FASTP_R2}" | \
    samtools view --threads "${SAMTOOLS_THREADS}" -bhS -o "${BWA_BAM}" -

  # sort and index
  samtools sort --threads "${SAMTOOLS_THREADS}" -o "${SAMTOOLS_SORTED}" -T "${BWA_ALIGN_DIR}/${SAMPLE}" "${BWA_BAM}"
  samtools index "${SAMTOOLS_SORTED}"
  # optionally remove unsorted bam
  rm -f "${BWA_BAM}"

  # Step 7: Trim alignments from amplicon scheme using align_trim.py
  echolog "Step 7: trim alignments (amplicon scheme) for ${SAMPLE}"
  mkdir -p "${RESULTS_DIR}/primertrimmed"
  PRIMER_UNSORTED_BAM="${RESULTS_DIR}/primertrimmed/${SAMPLE}.primertrimmed.unsorted.bam"
  PRIMER_SORTED_BAM="${RESULTS_DIR}/primertrimmed/${SAMPLE}.primertrimmed.rg.sorted.bam"
  ALIGNREPORT="${RESULTS_DIR}/primertrimmed/${SAMPLE}.alignreport.csv"
  AMP_DEPTH="${RESULTS_DIR}/primertrimmed/${SAMPLE}.amplicon_depths.tsv"
  ALIGN_ERR="${RESULTS_DIR}/primertrimmed/${SAMPLE}.alignreport.er"

  # run align_trim.py reading the sorted BAM as stdin and writing an unsorted BAM file
  # Step 7.1: Run align_trim.py to trim primers and normalize reads
  # - Takes a BAM file as input
  # - Uses the primer scheme (BED file)
  # - Outputs trimmed alignments in SAM format to a temporary file
  # - Writes alignment report and amplicon depth report
  python ./scripts/align_trim.py \
    --normalise "${NORMALISE}" \
    "${PRIMER_SCHEME_OUT}/primer.bed" \
    --paired \
    --no-read-groups \
    --primer-match-threshold "${PRIMER_MATCH_THRESHOLD}" \
    --min-mapq "${MIN_MAPQ}" \
    --trim-primers \
    --report "${ALIGNREPORT}" \
    --amp-depth-report "${AMP_DEPTH}" \
    < "${SAMTOOLS_SORTED}" \
    > "${PRIMER_UNSORTED_BAM}" \
    2> "${ALIGN_ERR}"

  # Sort + index the primer-trimmed BAM
  # Step 7.2: Convert and sort the SAM file into a BAM file
  samtools sort -T "${SAMPLE}" -o "${PRIMER_SORTED_BAM}" "${PRIMER_UNSORTED_BAM}"
  # Step 7.3: Index the sorted BAM file for fast random access
  samtools index "${PRIMER_SORTED_BAM}"
  rm -f "${PRIMER_UNSORTED_BAM}"

  # Step 8: Call variants with freebayes -> gVCF
  echolog "Step 8: calling variants (freebayes) for ${SAMPLE}"
  mkdir -p "${RESULTS_DIR}/freebayes"
  GVCF="${RESULTS_DIR}/freebayes/${SAMPLE}.gvcf"
  FREEBAYES_VARS="${RESULTS_DIR}/freebayes/${SAMPLE}.variants.vcf"
  FREEBAYES_CONS="${RESULTS_DIR}/freebayes/${SAMPLE}.consensus.vcf"
  FREEBAYES_MASK="${RESULTS_DIR}/freebayes/${SAMPLE}.mask.txt"

  freebayes \
    -p 1 \
    -f "${PRIMER_SCHEME_OUT}/reference.fasta" \
    -F "${FREEBAYES_FREQ_LOW}" \
    -C "${FREEBAYES_MIN_COUNT}" \
    --pooled-continuous \
    --min-coverage "${FREEBAYES_MIN_COVERAGE}" \
    --gvcf \
    --gvcf-dont-use-chunk true \
    "${PRIMER_SORTED_BAM}" \
    > "${GVCF}"

  # process gvcf -> mask, variants, consensus (script provided in tutorial)
  python ./scripts/process_gvcf.py \
    -d 10 \
    -l 0.25 \
    -u 0.75 \
    -m "${FREEBAYES_MASK}" \
    -v "${FREEBAYES_VARS}" \
    -c "${FREEBAYES_CONS}" \
    "${GVCF}"

  # Normalize both variant files using bcftools
  echolog "Normalizing VCFs for ${SAMPLE}"
  for v in "variants" "consensus"; do
    in_vcf="${RESULTS_DIR}/freebayes/${SAMPLE}.${v}.vcf"
    out_vcf="${RESULTS_DIR}/freebayes/${SAMPLE}.${v}.norm.vcf"
    bcftools norm -f "${PRIMER_SCHEME_OUT}/reference.fasta" "${in_vcf}" -O v -o "${out_vcf}"
    # bgzip and tabix the normalized file for later steps
    bgzip -f "${out_vcf}"
    tabix -f -p vcf "${out_vcf}.gz"
  done

  # Step 9: Generate consensus genome (ambiguous then fixed)
  echolog "Step 9: generating consensus for ${SAMPLE}"
  # split consensus.norm.vcf into ambiguous (IUPAC) and fixed
  CONS_NORM_GZ="${RESULTS_DIR}/freebayes/${SAMPLE}.consensus.norm.vcf.gz"

  for vt in "ambiguous" "fixed"; do
    out="${RESULTS_DIR}/freebayes/${SAMPLE}.${vt}.norm.vcf"
    gz_out="${out}.gz"
    awk -v vartag="ConsensusTag=${vt}" '$0 ~ /^#/ || $0 ~ vartag' <(zcat "${CONS_NORM_GZ}") > "${out}"
    bgzip -f "${out}"
    tabix -f -p vcf "${gz_out}"
  done

  # Apply ambiguous variants (IUPAC) first
  AMBIG_FA="${RESULTS_DIR}/freebayes/${SAMPLE}.ambiguous.fa"
  bcftools consensus -f "${PRIMER_SCHEME_OUT}/reference.fasta" -I "${RESULTS_DIR}/freebayes/${SAMPLE}.ambiguous.norm.vcf.gz" > "${AMBIG_FA}"

  # get contig name in reference (header without >)
  CTG_NAME=$(head -n1 "${PRIMER_SCHEME_OUT}/reference.fasta" | sed 's/>//')
  # produce final consensus applying fixed variants and mask, and rename contig header to sample name
  FINAL_CONS="${RESULTS_DIR}/freebayes/${SAMPLE}.consensus.fa"
  bcftools consensus -f "${AMBIG_FA}" -m "${FREEBAYES_MASK}" "${RESULTS_DIR}/freebayes/${SAMPLE}.fixed.norm.vcf.gz" | sed "s/${CTG_NAME}/${SAMPLE}/" > "${FINAL_CONS}"

  echolog "Sample ${SAMPLE} finished: consensus -> ${FINAL_CONS}"

  # Step 10: optionally run squirrel (user must set RUN_SQUIRREL=true at top)
  if [ "${RUN_SQUIRREL}" = "true" ]; then
    echolog "Step 10: running squirrel for ${SAMPLE} (may require additional sequences)"
    # Create/append consensus to all_consensus.fasta (build full set after processing all samples is better)
    cat "${FINAL_CONS}" >> "${RESULTS_DIR}/data/all-consensus/all_consensus.fasta"
  fi

  echolog "=== END SAMPLE: ${SAMPLE} ==="
}

### -----------------------
### Main loop: run for each sample
### -----------------------
for s in "${SAMPLES[@]}"; do
  process_sample "$s" || { echolog "ERROR processing sample $s"; exit 1; }
done

### -----------------------
### After sample loop: optional squirrel (Step 10 global)
### -----------------------
if [ "${RUN_SQUIRREL}" = "true" ]; then
  echolog "Running final Squirrel on concatenated consensus sequences"
  # ensure path to squirrel data & outdir exist
  mkdir -p "${RESULTS_DIR}/squirrel"
  # follow tutorial: assemble final upload file (requires additional fetch steps if you need them)
  # ensure file exists
  if [ ! -f "${RESULTS_DIR}/data/all-consensus/mpxv_all_consensus.fasta" ]; then
    # create by removing .1 suffixes if any
    sed 's/\.1.*//g' "${RESULTS_DIR}/data/all-consensus/all_consensus.fasta" > "${RESULTS_DIR}/data/all-consensus/mpxv_all_consensus.fasta"
  fi

  squirrel \
    "${RESULTS_DIR}/data/all-consensus/mpxv_all_consensus.fasta" \
    --no-mask \
    --seq-qc \
    --outdir "${RESULTS_DIR}/squirrel" \
    --outfile all_consensus.aln.fasta \
    --threads "${SQUIRREL_THREADS}" \
    --run-phylo \
    --run-apobec3-phylo \
    --interactive-tree \
    --outgroups KJ642617,KJ642615,KJ642616 \
    --clade cladei
fi

echolog "ALL SAMPLES PROCESSED. Script complete."

