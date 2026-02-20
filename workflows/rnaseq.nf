/*
========================================================================================
    RNASEQ-EDITING-QC: Main Workflow
========================================================================================
*/

nextflow.enable.dsl = 2

// Import nf-core modules
include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC           } from '../modules/nf-core/trimmomatic/main'
include { STAR_ALIGN            } from '../modules/nf-core/star/align/main'
include { HISAT2_ALIGN          } from '../modules/nf-core/hisat2/align/main'
include { SALMON_QUANT          } from '../modules/nf-core/salmon/quant/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'

// Import local custom modules
include { HTSEQ_COUNT           } from '../modules/local/htseq_count/main'
include { KO_VERIFICATION       } from '../modules/local/ko_verification/main'
include { DESEQ2_ANALYSIS       } from '../modules/local/deseq2_analysis/main'
include { GENERATE_REPORT       } from '../modules/local/generate_report/main'

// Import subworkflows
include { INPUT_CHECK           } from '../subworkflows/local/input_check'
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'

/*
========================================================================================
    WORKFLOW: RNASEQ_PIPELINE
========================================================================================
*/

workflow RNASEQ_PIPELINE {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Validate and parse samplesheet
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_reads = INPUT_CHECK.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Prepare genome indices
    //
    PREPARE_GENOME ()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // MODULE: FastQC on raw reads
    //
    if (!params.skip_fastqc) {
        FASTQC (ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //
    // MODULE: Trimming
    //
    ch_trimmed_reads = ch_reads
    if (!params.skip_trimming) {
        TRIMMOMATIC (ch_reads)
        ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
        ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    }

    //
    // MODULE: Alignment
    //
    ch_bam = Channel.empty()
    if (params.aligner == 'star') {
        STAR_ALIGN (
            ch_trimmed_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            false, // seq_platform
            false  // seq_center
        )
        ch_bam = STAR_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
    } else if (params.aligner == 'hisat2') {
        HISAT2_ALIGN (
            ch_trimmed_reads,
            PREPARE_GENOME.out.hisat2_index,
            PREPARE_GENOME.out.splicesites
        )
        ch_bam = HISAT2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
    }

    //
    // MODULE: Quantification
    //
    ch_counts = Channel.empty()
    if (params.quantification_method == 'htseq' || params.quantification_method == 'both') {
        HTSEQ_COUNT (
            ch_bam,
            PREPARE_GENOME.out.gtf
        )
        ch_counts = HTSEQ_COUNT.out.counts
        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
    }

    if (params.quantification_method == 'salmon' || params.quantification_method == 'both') {
        SALMON_QUANT (
            ch_trimmed_reads,
            PREPARE_GENOME.out.salmon_index,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.transcript_fasta,
            false, // alignment_mode
            params.salmon_quant_libtype
        )
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
    }

    //
    // MODULE: Knockout verification
    //
    KO_VERIFICATION (
        ch_counts.collect(),
        params.ko_gene,
        params.ko_condition,
        params.control_condition,
        params.ko_efficiency_threshold
    )
    ch_multiqc_files = ch_multiqc_files.mix(KO_VERIFICATION.out.multiqc_data.ifEmpty([]))

    //
    // MODULE: DESeq2 differential expression
    //
    DESEQ2_ANALYSIS (
        ch_counts.collect(),
        file(params.input),
        params.ko_condition,
        params.control_condition,
        params.lfc_threshold,
        params.pval_threshold
    )
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_ANALYSIS.out.multiqc_data.ifEmpty([]))
    ch_versions = ch_versions.mix(DESEQ2_ANALYSIS.out.versions)

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        MULTIQC (
            ch_multiqc_files.collect(),
            params.multiqc_config ? file(params.multiqc_config) : [],
            [],
            []
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    //
    // MODULE: Generate HTML report
    //
    if (!params.skip_report) {
        GENERATE_REPORT (
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
