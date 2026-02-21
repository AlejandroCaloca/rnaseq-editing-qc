#!/bin/bash
# =============================================================================
# install_nfcore_modules.sh
# Run this from INSIDE your cloned rnaseq-editing-qc repo directory.
#
# Prerequisites:
#   pipx install nf-core        (if not already installed)
#   git must be initialized    (already done)
#
# Usage:
#   cd ~/path/to/rnaseq-editing-qc
#   chmod +x install_nfcore_modules.sh
#   ./install_nfcore_modules.sh
# =============================================================================

set -e

echo "=============================================="
echo "  Installing nf-core modules"
echo "=============================================="

# Check nf-core is available
if ! command -v nf-core &> /dev/null; then
    echo "nf-core not found. Installing..."
    pipx install nf-core
    pipx ensurepath
fi

echo ""
echo "► QC modules"
nf-core modules install fastqc
nf-core modules install multiqc

echo ""
echo "► Alignment modules"
nf-core modules install star/genomegenerate
nf-core modules install star/align
nf-core modules install hisat2/build
nf-core modules install hisat2/align
nf-core modules install hisat2/extractsplicesites

echo ""
echo "► Quantification modules"
nf-core modules install salmon/index
nf-core modules install salmon/quant

echo ""
echo "► Trimming"
nf-core modules install trimmomatic

echo ""
echo "► BAM QC (samtools suite)"
nf-core modules install samtools/flagstat
nf-core modules install samtools/idxstats
nf-core modules install samtools/stats
nf-core modules install samtools/index
nf-core modules install samtools/sort

echo ""
echo "► RSeQC (strandedness inference)"
nf-core modules install rseqc/inferexperiment
nf-core modules install rseqc/readdistribution

echo ""
echo "=============================================="
echo "  All modules installed successfully!"
echo ""
echo "  Next steps:"
echo "  1. git add modules/"
echo "  2. git commit -m 'feat: install nf-core modules'"
echo "  3. git push"
echo "=============================================="
