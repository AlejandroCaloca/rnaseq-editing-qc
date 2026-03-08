/*
========================================================================================
    ENTRY POINT: BAM
    Pipeline starting from pre-aligned BAM files.
    Steps: BAM QC → Quantification → KO Verification → DESeq2 → Report
========================================================================================
*/

nextflow.enable.dsl = 2

include { INPUT_CHECK_BAM   } from '../subworkflows/local/input_check_bam'
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { BAM_QC            } from '../subworkflows/local/bam_qc'
include { SAMTOOLS_SORT     } from '../modules/nf-core/samtools/sort/main'
include { HTSEQ_COUNT       } from '../modules/local/htseq_count/main'
include { SALMON_QUANT      } from '../modules/nf-core/salmon/quant/main'
include { KO_VERIFICATION   } from '../modules/local/ko_verification/main'
include { DESEQ2_ANALYSIS   } from '../modules/local/deseq2_analysis/main'
include { GENERATE_REPORT   } from '../modules/local/generate_report/main'
include { MULTIQC           } from '../modules/nf-core/multiqc/main'

workflow PIPELINE_FROM_BAM {

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Validate BAM samplesheet → [ meta, bam, bai ]
    //
    INPUT_CHECK_BAM(file(params.input))
    ch_bam_input = INPUT_CHECK_BAM.out.bam

    //
    // SUBWORKFLOW: Prepare genome (GTF only for BAM entry point)
    //
    PREPARE_GENOME()

    //
    // SUBWORKFLOW: BAM QC
    //
    if (!params.skip_bam_qc) {
        BAM_QC(ch_bam_input, PREPARE_GENOME.out.gene_bed)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.multiqc_files.ifEmpty([]))
    }

    // Reshape for HTSeq: drop bai → [ meta, bam ]
    ch_bam = ch_bam_input.map { meta, bam, bai -> [ meta, bam ] }

    //
    // MODULE: Quantification
    //
    ch_counts = Channel.empty()

    if (params.quantification_method in ['htseq', 'both']) {
        // 1) Name-sort BAMS for HTseq --- HTSeq 2.0.2 can handle unsorted BAMs, but name-sorting is still recommended for better performance
        
        ch_bam_for_htseq = ch_bam

         // dummy optional inputs required by SAMTOOLS_SORT signature
         ch_no_fasta = ch_bam_for_htseq.map { meta,bam -> [meta, [null]] }   // tuple(meta2, fasta) with empty path
         ch_no_index = Channel.value(null)                      // val index_format

        SAMTOOLS_SORT(
            ch_bam_for_htseq,             // tuple [ meta, bam ]
            ch_no_fasta,       // tuple [ meta2, fasta ] with empty path
            ch_no_index        // val index_format with null value
                       
        )

        // 2) Map SAMTOOLS_SORT output to HTSEQ_COUNT input: [ meta, namesorted_bam ]
        ch_bam_namesorted = SAMTOOLS_SORT.out.bam.map { meta, bam_sorted -> [ meta, bam_sorted ] }

        // 3) Run HTSeq-count
        HTSEQ_COUNT(ch_bam, PREPARE_GENOME.out.gtf)
        ch_counts   = HTSEQ_COUNT.out.counts
    }

    if (params.quantification_method in ['salmon', 'both']) {
        SALMON_QUANT(
            ch_bam,
            PREPARE_GENOME.out.salmon_index,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.transcript_fasta,
            true,
            params.salmon_quant_libtype
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            SALMON_QUANT.out.results.collect{ it[1] }.ifEmpty([])
        )
    }

    //
    // MODULE: KO verification
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
    // MODULE: DESeq2
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

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        MULTIQC(
            ch_multiqc_files.collect(),
            params.multiqc_config ? file(params.multiqc_config) : [],
            [], []
        )
    }

    //
    // MODULE: HTML report
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
