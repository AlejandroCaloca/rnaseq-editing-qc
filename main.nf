#!/usr/bin/env nextflow

/*
========================================================================================
    rnaseq-editing-qc
========================================================================================
    IGF2BP3 Knockout RNA-seq Analysis Pipeline
    Custom fork of nf-core/rnaseq
    Author: Bioinformatics Research Lab
========================================================================================
*/

nextflow.enable.dsl = 2

// ── Help message ─────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
========================================================================================
    R N A S E Q - E D I T I N G - Q C   P I P E L I N E
========================================================================================

USAGE:
    nextflow run main.nf -profile <docker|singularity> --input <samplesheet.csv> [options]

ENTRY POINTS  (--entry_point):
    fastq          Start from raw FASTQ files  [DEFAULT]
                   Input columns: sample, fastq_1, fastq_2, condition, replicate
                   Steps: FastQC → Trimmomatic → STAR/HISAT2 → Quantification → DESeq2

    bam            Start from pre-aligned BAM files
                   Input columns: sample, bam, bai, condition, replicate
                   Steps: BAM QC (samtools/RSeQC) → Quantification → DESeq2

    counts         Start from a pre-computed count matrix or per-sample count files
                   Input columns: sample, counts_file, condition, replicate
                   Steps: DESeq2 + visualizations + report

    deseq2         Re-run plots and report from a saved DESeq2 RDS object
                   Use: --deseq2_rds /path/to/dds.rds
                   Steps: Plots + HTML report only

MANDATORY ARGS:
    --input        Path to samplesheet CSV (format depends on --entry_point)
    --genome       iGenomes key (e.g. GRCh38) OR provide --fasta + --gtf
    --outdir       Output directory [default: ./results]

COMMON OPTIONS:
    --entry_point       fastq | bam | counts | deseq2  [default: fastq]
    --ko_gene           Gene knocked out [default: IGF2BP3]
    --ko_condition      KO condition label in samplesheet [default: I3KO]
    --control_condition Control label [default: NT]
    --aligner           star | hisat2 [default: star]
    --lfc_threshold     |log2FC| cutoff [default: 1.0]
    --pval_threshold    padj cutoff [default: 0.05]
    --help              Show this help and exit

EXAMPLES:
    # Run from BAM files (most common current use case):
    nextflow run main.nf -profile docker \\
        --entry_point bam \\
        --input bam_samplesheet.csv \\
        --genome GRCh38 \\
        --outdir ./results

    # Run from raw FASTQs (full pipeline):
    nextflow run main.nf -profile docker \\
        --entry_point fastq \\
        --input fastq_samplesheet.csv \\
        --genome GRCh38 \\
        --outdir ./results

    # Re-run DE analysis from existing count files:
    nextflow run main.nf -profile docker \\
        --entry_point counts \\
        --input counts_samplesheet.csv \\
        --outdir ./results_rerun

    # Re-render plots/report only from saved DESeq2 RDS:
    nextflow run main.nf -profile docker \\
        --entry_point deseq2 \\
        --deseq2_rds ./results/deseq2/deseq2_object.rds \\
        --input samplesheet.csv \\
        --outdir ./results_replot

========================================================================================
""".stripIndent()
    exit 0
}

// ── Validate entry point ──────────────────────────────────────────────────────
def valid_entry_points = ['fastq', 'bam', 'counts', 'deseq2']
if (!valid_entry_points.contains(params.entry_point)) {
    exit 1, "ERROR: Invalid --entry_point '${params.entry_point}'.\n" +
             "       Valid options: ${valid_entry_points.join(', ')}\n" +
             "       Run with --help for usage examples."
}

// ── Validate inputs per entry point ──────────────────────────────────────────
if (params.entry_point == 'deseq2') {
    if (!params.deseq2_rds) {
        exit 1, "ERROR: --entry_point deseq2 requires --deseq2_rds /path/to/dds.rds"
    }
} else {
    if (!params.input) {
        exit 1, "ERROR: --input samplesheet not specified. Run with --help for usage."
    }
}

if (params.entry_point in ['fastq', 'bam'] && !params.genome && !params.fasta) {
    exit 1, "ERROR: --entry_point '${params.entry_point}' requires --genome or --fasta + --gtf"
}

// ── Print run header ──────────────────────────────────────────────────────────
log.info """
========================================================================================
    R N A S E Q - E D I T I N G - Q C
========================================================================================
    Entry point       : ${params.entry_point.toUpperCase()}
    Input             : ${params.input ?: 'N/A (deseq2 RDS mode)'}
    Genome            : ${params.genome ?: params.fasta ?: 'N/A'}
    Output dir        : ${params.outdir}
    KO gene           : ${params.ko_gene}
    Contrast          : ${params.ko_condition} vs ${params.control_condition}
    Aligner           : ${params.aligner}
    Quantification    : ${params.quantification_method}
    LFC threshold     : ${params.lfc_threshold}
    Padj threshold    : ${params.pval_threshold}
========================================================================================
""".stripIndent()

// ── Import entry-point workflows ──────────────────────────────────────────────
include { PIPELINE_FROM_FASTQ  } from './workflows/entry_fastq'
include { PIPELINE_FROM_BAM    } from './workflows/entry_bam'
include { PIPELINE_FROM_COUNTS } from './workflows/entry_counts'
include { PIPELINE_FROM_DESEQ2 } from './workflows/entry_deseq2'

/*
========================================================================================
    MAIN — routes execution to the correct entry-point workflow
========================================================================================
*/

workflow {
    switch (params.entry_point) {
        case 'fastq':
            PIPELINE_FROM_FASTQ()
            break
        case 'bam':
            PIPELINE_FROM_BAM()
            break
        case 'counts':
            PIPELINE_FROM_COUNTS()
            break
        case 'deseq2':
            PIPELINE_FROM_DESEQ2()
            break
    }
}

workflow.onComplete {
    def status = workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'
    log.info """
========================================================================================
    Pipeline ${status}
    Entry point : ${params.entry_point.toUpperCase()}
    Results     : ${params.outdir}
    Duration    : ${workflow.duration}
========================================================================================
""".stripIndent()
}
