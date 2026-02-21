#!/bin/bash
# =============================================================================
# index_bams_and_test.sh
#
# Step 1: Indexes all BAM files in a directory using samtools
# Step 2: Auto-generates a BAM samplesheet from the indexed BAMs
# Step 3: Runs the pipeline with --entry_point bam
#
# Usage:
#   chmod +x index_bams_and_test.sh
#   ./index_bams_and_test.sh /path/to/bam/directory
#
# Requirements:
#   - samtools (install: sudo apt install samtools  OR  conda install samtools)
#   - nextflow installed
#   - docker running (or use -profile singularity)
# =============================================================================

set -e

# ── Arguments ─────────────────────────────────────────────────────────────────
BAM_DIR="${1:-.}"          # Directory containing BAM files (default: current dir)
OUTDIR="./results_bam_test"
GENOME="GRCh38"            # Change to your genome key, or set FASTA/GTF below
PROFILE="docker"           # Change to 'singularity' if on HPC

# If using a custom genome instead of iGenomes, set these:
FASTA=""   # e.g. /data/genome/GRCh38.fa
GTF=""     # e.g. /data/genome/GRCh38.gtf

# Pipeline parameters
KO_CONDITION="I3KO"
CONTROL_CONDITION="NT"
KO_GENE="IGF2BP3"

echo "=============================================="
echo "  BAM Indexing + Pipeline Test"
echo "  BAM directory: ${BAM_DIR}"
echo "=============================================="

# ── Step 1: Check samtools ────────────────────────────────────────────────────
if ! command -v samtools &> /dev/null; then
    echo ""
    echo "ERROR: samtools not found."
    echo "Install with one of:"
    echo "  sudo apt install samtools"
    echo "  conda install -c bioconda samtools"
    exit 1
fi

SAMTOOLS_VERSION=$(samtools --version | head -1)
echo "Using: ${SAMTOOLS_VERSION}"

# ── Step 2: Find and index BAM files ─────────────────────────────────────────
echo ""
echo "► Scanning for BAM files in: ${BAM_DIR}"

BAM_FILES=($(find "${BAM_DIR}" -name "*.bam" | sort))

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "ERROR: No .bam files found in ${BAM_DIR}"
    echo "Usage: ./index_bams_and_test.sh /path/to/bam/folder"
    exit 1
fi

echo "  Found ${#BAM_FILES[@]} BAM file(s):"
for bam in "${BAM_FILES[@]}"; do
    echo "    - $(basename ${bam})"
done

echo ""
echo "► Indexing BAM files..."
for bam in "${BAM_FILES[@]}"; do
    bai="${bam}.bai"
    if [ -f "${bai}" ]; then
        echo "  [SKIP] Index already exists: $(basename ${bai})"
    else
        echo "  [INDEX] $(basename ${bam}) ..."
        samtools index -@ 4 "${bam}"
        echo "  [DONE]  $(basename ${bai}) created"
    fi
done

echo ""
echo "► All BAM files indexed."

# ── Step 3: Auto-generate samplesheet ────────────────────────────────────────
echo ""
echo "► Generating BAM samplesheet..."

SAMPLESHEET="./bam_samplesheet_auto.csv"
echo "sample,bam,bai,condition,replicate" > "${SAMPLESHEET}"

REP_KO=1
REP_NT=1

for bam in "${BAM_FILES[@]}"; do
    bam_abs=$(realpath "${bam}")
    bai_abs="${bam_abs}.bai"
    filename=$(basename "${bam}" .bam)

    # Infer condition from filename
    # Looks for KO_CONDITION or CONTROL_CONDITION string in filename
    if echo "${filename}" | grep -qi "${KO_CONDITION}"; then
        condition="${KO_CONDITION}"
        replicate="${REP_KO}"
        REP_KO=$((REP_KO + 1))
        sample_name="SEM_${KO_CONDITION}_rep${replicate}"
    elif echo "${filename}" | grep -qi "${CONTROL_CONDITION}"; then
        condition="${CONTROL_CONDITION}"
        replicate="${REP_NT}"
        REP_NT=$((REP_NT + 1))
        sample_name="SEM_${CONTROL_CONDITION}_rep${replicate}"
    else
        echo "  [WARN] Could not infer condition from filename: ${filename}"
        echo "         Filename should contain '${KO_CONDITION}' or '${CONTROL_CONDITION}'"
        echo "         Defaulting to UNKNOWN — edit the samplesheet manually."
        condition="UNKNOWN"
        replicate=1
        sample_name="${filename}"
    fi

    echo "${sample_name},${bam_abs},${bai_abs},${condition},${replicate}" >> "${SAMPLESHEET}"
    echo "  Added: ${sample_name} → ${condition} (rep ${replicate})"
done

echo ""
echo "  Samplesheet written to: ${SAMPLESHEET}"
echo ""
echo "  ── Contents ──────────────────────────────"
cat "${SAMPLESHEET}"
echo "  ──────────────────────────────────────────"
echo ""
echo "  IMPORTANT: Review the samplesheet before running!"
echo "  If conditions or replicates look wrong, edit ${SAMPLESHEET} manually."
echo ""
read -p "  Does the samplesheet look correct? (y/n): " confirm
if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
    echo "  Exiting. Edit ${SAMPLESHEET} and re-run the pipeline manually (see below)."
    echo ""
    print_manual_run_cmd
    exit 0
fi

# ── Step 4: Run pipeline ──────────────────────────────────────────────────────
echo ""
echo "► Launching pipeline..."
echo ""

# Build genome argument
if [ -n "${FASTA}" ] && [ -n "${GTF}" ]; then
    GENOME_ARGS="--fasta ${FASTA} --gtf ${GTF}"
else
    GENOME_ARGS="--genome ${GENOME}"
fi

nextflow run main.nf \
    -profile "${PROFILE}" \
    --entry_point bam \
    --input "${SAMPLESHEET}" \
    ${GENOME_ARGS} \
    --ko_gene "${KO_GENE}" \
    --ko_condition "${KO_CONDITION}" \
    --control_condition "${CONTROL_CONDITION}" \
    --outdir "${OUTDIR}" \
    -resume

echo ""
echo "=============================================="
echo "  Pipeline complete!"
echo "  Results: ${OUTDIR}"
echo "  Report:  ${OUTDIR}/report/*.html"
echo "=============================================="

# Helper function
print_manual_run_cmd() {
    echo "  Manual run command:"
    echo ""
    echo "    nextflow run main.nf \\"
    echo "        -profile ${PROFILE} \\"
    echo "        --entry_point bam \\"
    echo "        --input ${SAMPLESHEET} \\"
    echo "        --genome ${GENOME} \\"
    echo "        --ko_gene ${KO_GENE} \\"
    echo "        --ko_condition ${KO_CONDITION} \\"
    echo "        --control_condition ${CONTROL_CONDITION} \\"
    echo "        --outdir ${OUTDIR}"
    echo ""
}
