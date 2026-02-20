/*
========================================================================================
    ENTRY POINT: BAM
    Pipeline starting from pre-aligned BAM files — skips FastQC, trimming, alignment.
    Steps: BAM QC → Quantification → KO Verification → DESeq2 → Report

    SAMPLESHEET FORMAT (bam_samplesheet.csv):
        sample,bam,bai,condition,replicate
        SEM_I3KO_rep1,/path/to/I3KO_rep1.bam,/path/to/I3KO_rep1.bam.bai,I3KO,1
        SEM_NT_rep1,/path/to/NT_rep1.bam,/path/to/NT_rep1.bam.bai,NT,1
        ...
========================================================================================
*/

nextflow.enable.dsl = 2

include { INPUT_CHECK_BAM       } from '../subworkflows/local/input_check_bam'
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'
include { BAM_QC                } from '../subworkflows/local/bam_qc'
include { HTSEQ_COUNT           } from '../modules/local/htseq_count/main'
include { SALMON_QUANT          } from '../modules/nf-core/salmon/quant/main'
include { KO_VERIFICATION       } from '../modules/local/ko_verification/main'
include { DESEQ2_ANALYSIS       } from '../modules/local/deseq2_analysis/main'
include { GENERATE_REPORT       } from '../modules/local/generate_report/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'

workflow PIPELINE_FROM_BAM {

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Validate BAM samplesheet → emit [ meta, bam, bai ] tuples
    //
    INPUT_CHECK_BAM(file(params.input))
    ch_bam_input = INPUT_CHECK_BAM.out.bam   // [ meta, bam, bai ]
    ch_versions  = ch_versions.mix(INPUT_CHECK_BAM.out.versions)

    //
    // SUBWORKFLOW: Prepare genome (GTF needed for quantification)
    //
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: BAM QC — samtools flagstat/idxstats + RSeQC (optional)
    //
    if (!params.skip_bam_qc) {
        BAM_QC(ch_bam_input, PREPARE_GENOME.out.gene_bed)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.multiqc_files.ifEmpty([]))
        ch_versions      = ch_versions.mix(BAM_QC.out.versions)
    }

    // Reshape for downstream: [ meta, bam ]  (drop bai — not needed by htseq)
    ch_bam = ch_bam_input.map { meta, bam, bai -> [ meta, bam ] }

    //
    // MODULE: Quantification
    //
    ch_counts = Channel.empty()

    if (params.quantification_method in ['htseq', 'both']) {
        HTSEQ_COUNT(ch_bam, PREPARE_GENOME.out.gtf)
        ch_counts   = HTSEQ_COUNT.out.counts
        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
    }

    if (params.quantification_method in ['salmon', 'both']) {
        // Salmon alignment-based mode: accepts BAM directly
        SALMON_QUANT(
            ch_bam,
            PREPARE_GENOME.out.salmon_index,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.transcript_fasta,
            true,   // alignment_mode = true
            params.salmon_quant_libtype
        )
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{ it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(SALMON_QUANT.out.versions.first())
    }

    //
    // MODULE: KO verification (checks IGF2BP3 is knocked out)
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
    // MODULE: DESeq2 differential expression
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
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        MULTIQC(
            ch_multiqc_files.collect(),
            params.multiqc_config ? file(params.multiqc_config) : [],
            [], []
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    //
    // MODULE: Final HTML report
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
