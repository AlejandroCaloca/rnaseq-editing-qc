/*
========================================================================================
    ENTRY POINT: FASTQ
    Full pipeline starting from raw paired-end FASTQ files
    Steps: FastQC → Trimming → Alignment → Quantification → KO Check → DESeq2 → Report
========================================================================================
*/

nextflow.enable.dsl = 2

include { INPUT_CHECK_FASTQ     } from '../subworkflows/local/input_check_fastq'
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'
include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC           } from '../modules/nf-core/trimmomatic/main'
include { STAR_ALIGN            } from '../modules/nf-core/star/align/main'
include { HISAT2_ALIGN          } from '../modules/nf-core/hisat2/align/main'
include { HTSEQ_COUNT           } from '../modules/local/htseq_count/main'
include { SALMON_QUANT          } from '../modules/nf-core/salmon/quant/main'
include { KO_VERIFICATION       } from '../modules/local/ko_verification/main'
include { DESEQ2_ANALYSIS       } from '../modules/local/deseq2_analysis/main'
include { GENERATE_REPORT       } from '../modules/local/generate_report/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'

workflow PIPELINE_FROM_FASTQ {

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Validate samplesheet and create reads channel
    INPUT_CHECK_FASTQ(file(params.input))
    ch_reads    = INPUT_CHECK_FASTQ.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK_FASTQ.out.versions)

    // Prepare genome (build or reuse indices)
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Raw QC
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(FASTQC.out.versions.first())
    }

    // Trimming
    ch_trimmed = ch_reads
    if (!params.skip_trimming) {
        TRIMMOMATIC(ch_reads)
        ch_trimmed       = TRIMMOMATIC.out.trimmed_reads
        ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{ it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    }

    // Alignment
    ch_bam = Channel.empty()
    if (params.aligner == 'star') {
        STAR_ALIGN(ch_trimmed, PREPARE_GENOME.out.star_index, PREPARE_GENOME.out.gtf, false, false)
        ch_bam           = STAR_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{ it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(STAR_ALIGN.out.versions.first())
    } else {
        HISAT2_ALIGN(ch_trimmed, PREPARE_GENOME.out.hisat2_index, PREPARE_GENOME.out.splicesites)
        ch_bam           = HISAT2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{ it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
    }

    // Route to shared quantification + downstream
    _run_downstream(ch_bam, ch_multiqc_files, ch_versions, PREPARE_GENOME.out.gtf,
                    PREPARE_GENOME.out.salmon_index, PREPARE_GENOME.out.transcript_fasta)
}

// Shared block reused by both fastq and bam entry points
def _run_downstream(ch_bam, ch_multiqc_files, ch_versions, ch_gtf, ch_salmon_index, ch_transcript_fasta) {

    // Quantification
    ch_counts = Channel.empty()
    if (params.quantification_method in ['htseq', 'both']) {
        HTSEQ_COUNT(ch_bam, ch_gtf)
        ch_counts   = HTSEQ_COUNT.out.counts
        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
    }
    if (params.quantification_method in ['salmon', 'both']) {
        // Salmon in alignment-based mode using BAM
        SALMON_QUANT(ch_bam, ch_salmon_index, ch_gtf, ch_transcript_fasta, true, params.salmon_quant_libtype)
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{ it[1] }.ifEmpty([]))
        ch_versions      = ch_versions.mix(SALMON_QUANT.out.versions.first())
    }

    // KO verification
    KO_VERIFICATION(
        ch_counts.collect(),
        params.ko_gene, params.ko_condition, params.control_condition,
        params.ko_efficiency_threshold
    )
    ch_multiqc_files = ch_multiqc_files.mix(KO_VERIFICATION.out.multiqc_data.ifEmpty([]))

    // Differential expression
    DESEQ2_ANALYSIS(
        ch_counts.collect(), file(params.input),
        params.ko_condition, params.control_condition,
        params.lfc_threshold, params.pval_threshold
    )
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_ANALYSIS.out.multiqc_data.ifEmpty([]))
    ch_versions      = ch_versions.mix(DESEQ2_ANALYSIS.out.versions)

    // MultiQC
    if (!params.skip_multiqc) {
        MULTIQC(ch_multiqc_files.collect(), params.multiqc_config ? file(params.multiqc_config) : [], [], [])
    }

    // Final report
    if (!params.skip_report) {
        GENERATE_REPORT(
            DESEQ2_ANALYSIS.out.results_all, DESEQ2_ANALYSIS.out.results_sig,
            KO_VERIFICATION.out.ko_stats, KO_VERIFICATION.out.ko_plot,
            DESEQ2_ANALYSIS.out.plots, file(params.input), ch_versions.collect()
        )
    }
}
