#!/usr/bin/env bash
set -euo pipefail

# -------- CONFIG --------
BAM_IN="sem_i3ko_dmso_zn264_hg38.ih1.uddD1.xpe.yc.star.map.trk.bam"
GTF="gencode.v47.primary_assembly.annotation.gtf"
OUTDIR="./results_test_6samples"
SAMPLE_CSV="assets/samplesheet_bam_test_6samples.csv"

# Subsample fractions for pseudo-replicates
SUBSAMPLE_FRACS=(0.01 0.02 0.03 0.04 0.05 0.06)

# -------- MAKE 6 SMALL BAMs --------
echo "Subsampling input BAM into 6 small test BAMs..."
for i in "${!SUBSAMPLE_FRACS[@]}"; do
  idx=$((i+1))
  frac=${SUBSAMPLE_FRACS[$i]}
  out="test_rep${idx}.bam"
  samtools view -s ${frac} -b "${BAM_IN}" > "${out}"
  samtools index "${out}"
done

# -------- CREATE TEST SAMPLESHEET --------
cat > "${SAMPLE_CSV}" <<EOF
sample,bam,bai,condition,replicate
SEM_I3KO_rep1,test_rep1.bam,test_rep1.bam.bai,I3KO,1
SEM_I3KO_rep2,test_rep2.bam,test_rep2.bam.bai,I3KO,2
SEM_I3KO_rep3,test_rep3.bam,test_rep3.bam.bai,I3KO,3
SEM_NT_rep1,test_rep4.bam,test_rep4.bam.bai,NT,1
SEM_NT_rep2,test_rep5.bam,test_rep5.bam.bai,NT,2
SEM_NT_rep3,test_rep6.bam,test_rep6.bam.bai,NT,3
EOF

# -------- RUN NEXTFLOW --------
nextflow run main.nf \
  -profile docker \
  --entryPoint bam \
  --input "${SAMPLE_CSV}" \
  --gtf "${GTF}" \
  --skip_multiqc \
  --skip_report \
  --outdir "${OUTDIR}" \
  -resume