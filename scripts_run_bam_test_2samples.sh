#!/usr/bin/env bash
set -euo pipefail

# -------- CONFIG --------
BAM_EXP="path/to/experimental.bam"
BAM_CTRL="path/to/control.bam"
GTF="gencode.v49.primary_assembly.annotation.gtf"
OUTDIR="./results_test_2samples"
SAMPLE_CSV="assets/samplesheet_bam_test_2samples.csv"

# -------- VALIDATE --------
for f in "${BAM_EXP}" "${BAM_EXP}.bai" "${BAM_CTRL}" "${BAM_CTRL}.bai" "${GTF}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Missing file: ${f}"
    exit 1
  fi
done

# -------- CREATE SAMPLESHEET --------
cat > "${SAMPLE_CSV}" <<EOF
sample,bam,bai,condition,replicate
EXP_rep1,${BAM_EXP},${BAM_EXP}.bai,I3KO,1
CTRL_rep1,${BAM_CTRL},${BAM_CTRL}.bai,NT,1
EOF

# -------- RUN NEXTFLOW --------
nextflow run main.nf \
  -profile slurm,podman \
  --entryPoint bam \
  --input "${SAMPLE_CSV}" \
  --gtf "${GTF}" \
  --skip_multiqc \
  --skip_report \
  --outdir "${OUTDIR}" \
  -resume