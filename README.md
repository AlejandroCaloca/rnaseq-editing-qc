# rnaseq-editing-qc

> **RNA-seq pipeline for genome-edited cell line analysis**  
> IGF2BP3 knockout (I3KO) vs. Non-targeting control (NT) in SEM cells

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Summary](#pipeline-summary)
3. [Quick Start](#quick-start)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Samplesheet Format](#samplesheet-format)
7. [Parameters](#parameters)
8. [Outputs](#outputs)
9. [Configuration Profiles](#configuration-profiles)
10. [Running on HPC (SLURM)](#running-on-hpc-slurm)
11. [Troubleshooting](#troubleshooting)
12. [Citation](#citation)

---

## Overview

This pipeline is a custom fork of [nf-core/rnaseq](https://github.com/nf-core/rnaseq) designed specifically for **comparing CRISPR-edited knockout cell lines against non-targeting controls**. It automates the full analysis from raw FASTQ files to a publication-ready HTML report.

**Designed for:**
- SEM (B-cell precursor leukemia) cell line experiments
- IGF2BP3 knockout (I3KO) vs. non-targeting control (NT) comparisons
- 3 biological replicates per condition, paired-end sequencing

---

## Pipeline Summary

```
Raw FASTQs
    │
    ▼
FastQC (raw QC)
    │
    ▼
Trimmomatic (adapter & quality trimming)
    │
    ▼
FastQC (post-trim QC)
    │
    ▼
STAR or HISAT2 (alignment to GRCh38)
    │
    ▼
HTSeq-count / Salmon (quantification)
    │
    ├──► Knockout Verification (IGF2BP3 expression check)
    │
    ▼
DESeq2 (differential expression: I3KO vs NT)
    │
    ├──► Volcano plot
    ├──► MA plot
    ├──► PCA plot
    ├──► Heatmap (top DE genes)
    │
    ▼
MultiQC (aggregate QC report)
    │
    ▼
HTML Report (RMarkdown, main deliverable)
```

---

## Quick Start

```bash
# 1. Install Nextflow (requires Java 11+)
curl -s https://get.nextflow.io | bash

# 2. Clone this repository
git clone https://github.com/your-org/rnaseq-editing-qc.git
cd rnaseq-editing-qc

# 3. Run with Docker (recommended)
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --genome GRCh38 \
    --outdir ./results

# 4. Run test to validate installation
nextflow run main.nf -profile docker,test
```

---

## Installation

### Requirements

| Tool | Version | Notes |
|------|---------|-------|
| Java | ≥ 11 | Required for Nextflow |
| Nextflow | ≥ 22.10.0 | Pipeline manager |
| Docker **or** Singularity | Any | For containerized execution |

### Install Nextflow

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

### Clone Repository

```bash
git clone https://github.com/your-org/rnaseq-editing-qc.git
cd rnaseq-editing-qc
```

---

## Usage

### Basic Command

```bash
nextflow run main.nf \
    -profile docker \
    --input /path/to/samplesheet.csv \
    --genome GRCh38 \
    --outdir ./results
```

### Using a Custom Genome

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --fasta /path/to/GRCh38.fa \
    --gtf /path/to/GRCh38.gtf \
    --outdir ./results
```

### With Pre-built Indices (faster)

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --genome GRCh38 \
    --star_index /path/to/star_index/ \
    --outdir ./results
```

### Resuming a Failed Run

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --genome GRCh38 -resume
```

---

## Samplesheet Format

Create a CSV file with the following columns:

```
sample,fastq_1,fastq_2,condition,replicate
```

**Example** (`samplesheet.csv`):

```csv
sample,fastq_1,fastq_2,condition,replicate
SEM_I3KO_DMSO_rep1,/data/I3KO_rep1_R1.fastq.gz,/data/I3KO_rep1_R2.fastq.gz,I3KO,1
SEM_I3KO_DMSO_rep2,/data/I3KO_rep2_R1.fastq.gz,/data/I3KO_rep2_R2.fastq.gz,I3KO,2
SEM_I3KO_DMSO_rep3,/data/I3KO_rep3_R1.fastq.gz,/data/I3KO_rep3_R2.fastq.gz,I3KO,3
SEM_NT_DMSO_rep1,/data/NT_rep1_R1.fastq.gz,/data/NT_rep1_R2.fastq.gz,NT,1
SEM_NT_DMSO_rep2,/data/NT_rep2_R1.fastq.gz,/data/NT_rep2_R2.fastq.gz,NT,2
SEM_NT_DMSO_rep3,/data/NT_rep3_R1.fastq.gz,/data/NT_rep3_R2.fastq.gz,NT,3
```

A template is available at `assets/samplesheet_template.csv`.

**Column descriptions:**
- `sample`: Unique sample name (no spaces)
- `fastq_1`: Path to R1 FASTQ file (gzipped)
- `fastq_2`: Path to R2 FASTQ file (gzipped)
- `condition`: Condition label — must match `--ko_condition` or `--control_condition`
- `replicate`: Integer replicate number

---

## Parameters

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | **required** | Path to samplesheet CSV |
| `--genome` | null | iGenomes key (e.g. `GRCh38`) |
| `--fasta` | null | Path to genome FASTA (alternative to `--genome`) |
| `--gtf` | null | Path to GTF annotation file |
| `--outdir` | `./results` | Output directory |

### Experiment Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ko_gene` | `IGF2BP3` | Name of the knocked-out gene |
| `--ko_condition` | `I3KO` | Label for knockout condition in samplesheet |
| `--control_condition` | `NT` | Label for control condition in samplesheet |
| `--ko_efficiency_threshold` | `0.90` | Minimum KO efficiency to pass QC (0–1) |

### Alignment Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aligner` | `star` | Aligner to use: `star` or `hisat2` |
| `--star_index` | null | Pre-built STAR index directory |
| `--hisat2_index` | null | Pre-built HISAT2 index directory |

### DESeq2 Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--lfc_threshold` | `1.0` | `\|log2FC\|` cutoff for significance |
| `--pval_threshold` | `0.05` | Adjusted p-value cutoff |
| `--top_de_genes` | `50` | Number of top DE genes in heatmap |
| `--genes_of_interest` | `''` | Comma-separated genes for barplots |

### Skip Flags

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_fastqc` | false | Skip FastQC |
| `--skip_trimming` | false | Skip trimming |
| `--skip_multiqc` | false | Skip MultiQC |
| `--skip_report` | false | Skip HTML report generation |

---

## Outputs

```
results/
├── fastqc/                     # FastQC reports (raw + trimmed)
├── trimmomatic/                # Trimmed reads and logs
├── star/                       # Aligned BAM files
├── htseq/                      # Per-sample count files
├── salmon/                     # Salmon quantification
├── ko_verification/
│   ├── ko_efficiency_stats.csv # KO efficiency metrics
│   ├── ko_efficiency_plot.pdf  # Barplot of KO gene expression
│   └── ko_efficiency_plot.png
├── deseq2/
│   ├── deseq2_results_all.csv  # All genes, ranked by padj
│   ├── deseq2_results_all.xlsx # Same, Excel format
│   ├── deseq2_results_sig.csv  # Significant genes only
│   ├── deseq2_results_sig.xlsx
│   ├── deseq2_object.rds       # DESeq2 R object for further analysis
│   └── plots/
│       ├── pca_plot.pdf / .png
│       ├── volcano_plot.pdf / .png
│       ├── ma_plot.pdf / .png
│       └── heatmap_top_genes.pdf / .png
├── multiqc/
│   └── multiqc_report.html     # Aggregate QC report
├── report/
│   └── I3KO_vs_NT_analysis_report.html  # ⭐ Main deliverable
└── pipeline_info/              # Nextflow execution logs
```

---

## Configuration Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (recommended for local) |
| `singularity` | Run with Singularity (recommended for HPC) |
| `slurm` | Submit jobs via SLURM scheduler |
| `test` | Run with minimal test data to validate installation |
| `test_full` | Full-scale test run |

Profiles can be combined:
```bash
-profile singularity,slurm
```

---

## Running on HPC (SLURM)

```bash
# Submit pipeline as a SLURM job
sbatch << 'EOF'
#!/bin/bash
#SBATCH --job-name=rnaseq_igf2bp3
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=pipeline_%j.log

module load nextflow
module load singularity

nextflow run main.nf \
    -profile singularity,slurm \
    --input /data/samplesheet.csv \
    --genome GRCh38 \
    --outdir /data/results \
    --slurm_account YOUR_ACCOUNT \
    -resume
EOF
```

---

## Troubleshooting

**Pipeline fails at STAR alignment with memory error**
```bash
# Increase memory in nextflow.config or use --max_memory flag
nextflow run main.nf --max_memory 64.GB ...
```

**"Gene not found" error in KO verification**  
Check that `--ko_gene` exactly matches the gene name in your GTF (case-sensitive).

**Samples fail KO efficiency QC**  
The pipeline will warn but continue. Check the `ko_efficiency_stats.csv` file.  
Common causes: incomplete editing, contamination, or wrong condition labels in samplesheet.

**Resume not working**  
Delete the `.nextflow/` directory and `work/` folder, then rerun without `-resume`.

For further help, open an issue on GitHub.

---

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso et al., *Nature Biotechnology* (2017)
- **nf-core**: Ewels et al., *Nature Biotechnology* (2020)
- **DESeq2**: Love et al., *Genome Biology* (2014)
- **STAR**: Dobin et al., *Bioinformatics* (2013)
- **FastQC**: Andrews S. (2010) — https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: Ewels et al., *Bioinformatics* (2016)
- **HTSeq**: Putri et al., *Bioinformatics* (2022)
- **Trimmomatic**: Bolger et al., *Bioinformatics* (2014)

---

*Developed by the Bioinformatics Research Lab*
