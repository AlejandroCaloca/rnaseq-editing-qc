/*
========================================================================================
    ENTRY POINT: COUNTS
    Pipeline starting from pre-computed count files — skips all upstream processing.
    Steps: DESeq2 → KO Verification → Visualizations → Report

    SAMPLESHEET FORMAT — per-sample HTSeq files (counts_format = 'htseq'):
        sample,counts_file,condition,replicate
        SEM_I3KO_rep1,/path/to/I3KO_rep1.counts.txt,I3KO,1
        SEM_NT_rep1,/path/to/NT_rep1.counts.txt,NT,1

    SAMPLESHEET FORMAT — single count matrix (counts_format = 'matrix'):
        Provide --counts_matrix /path/to/counts_matrix.csv instead of --input
        Matrix format: rows=genes, cols=samples (first col = gene_id)
========================================================================================
*/

nextflow.enable.dsl = 2

include { INPUT_CHECK_COUNTS    } from '../subworkflows/local/input_check_counts'
include { KO_VERIFICATION       } from '../modules/local/ko_verification/main'
include { DESEQ2_ANALYSIS       } from '../modules/local/deseq2_analysis/main'
include { GENERATE_REPORT       } from '../modules/local/generate_report/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'

workflow PIPELINE_FROM_COUNTS {

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Validate samplesheet and emit count file paths
    //
    INPUT_CHECK_COUNTS(file(params.input))
    ch_counts   = INPUT_CHECK_COUNTS.out.counts   // [ meta, counts_file ]
    ch_versions = ch_versions.mix(INPUT_CHECK_COUNTS.out.versions)

    //
    // KO verification
    //
    KO_VERIFICATION(
        ch_counts.collect(),
        params.ko_gene,
        params.ko_condition,
        params.control_condition,
        params.ko_efficiency_threshold
    )
    ch_multiqc_files = ch_multiqc_files.mix(KO_VERIFICATION.out.multiqc_data.ifEmpty([]))

    //
    // DESeq2
    //
    DESEQ2_ANALYSIS(
        ch_counts.collect(),
        file(params.input),
        params.ko_condition,
        params.control_condition,
        params.lfc_threshold,
        params.pval_threshold
    )
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_ANALYSIS.out.multiqc_data.ifEmpty([]))
    ch_versions      = ch_versions.mix(DESEQ2_ANALYSIS.out.versions)

    //
    // MultiQC
    //
    if (!params.skip_multiqc) {
        MULTIQC(ch_multiqc_files.collect(), [], [], [])
    }

    //
    // Report
    //
    if (!params.skip_report) {
        GENERATE_REPORT(
            DESEQ2_ANALYSIS.out.results_all,
            DESEQ2_ANALYSIS.out.results_sig,
            KO_VERIFICATION.out.ko_stats,
            KO_VERIFICATION.out.ko_plot,
            DESEQ2_ANALYSIS.out.plots,
            file(params.input),
            ch_versions.collect()
        )
    }
}
