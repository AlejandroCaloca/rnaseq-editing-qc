#!/bin/bash
# =============================================================================
# setup_github.sh
# Run this script ONCE to initialize the local git repo and push to GitHub
# 
# Usage:
#   chmod +x setup_github.sh
#   ./setup_github.sh your-github-org-or-username
# =============================================================================

set -e

ORG="${1:-your-org}"
REPO="rnaseq-editing-qc"
REMOTE="https://github.com/${ORG}/${REPO}.git"

echo "==================================================="
echo "  Setting up GitHub repo: ${REMOTE}"
echo "==================================================="

# Initialize git repo
git init
git add .
git commit -m "feat: initial pipeline scaffold

- main.nf entry point with parameter validation
- nextflow.config with all tunable parameters
- workflows/rnaseq.nf - main workflow DAG
- subworkflows/local/input_check.nf - samplesheet validation
- subworkflows/local/prepare_genome.nf - index building
- modules/local/htseq_count/main.nf
- modules/local/ko_verification/main.nf - IGF2BP3 QC module
- modules/local/deseq2_analysis/main.nf - DE analysis + plots
- modules/local/generate_report/main.nf - HTML report
- conf/base.config, test.config, slurm.config
- assets/samplesheet_template.csv
- test_data/samplesheet_test.csv
- README.md with full documentation
- .gitignore"

git branch -M main

echo ""
echo "Next steps:"
echo "  1. Create the repo on GitHub: https://github.com/new"
echo "     Name: ${REPO}"
echo "     Visibility: Private (recommended)"
echo "     Do NOT initialize with README (we already have one)"
echo ""
echo "  2. Then run:"
echo "     git remote add origin ${REMOTE}"
echo "     git push -u origin main"
echo ""
echo "  3. Set up branch protection on 'main' in GitHub Settings > Branches"
echo ""
echo "Suggested GitHub repo labels to create:"
echo "  - enhancement (blue)"
echo "  - bug (red)"
echo "  - documentation (green)"
echo "  - module: qc"
echo "  - module: deseq2"
echo "  - module: report"
echo "  - priority: high"
