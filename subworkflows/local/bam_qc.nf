/*
========================================================================================
    SUBWORKFLOW: BAM_QC
    Quality control of pre-aligned BAM files.
    Runs: samtools flagstat, samtools idxstats
    Note: samtools stats requires a FASTA reference and is excluded for BAM entry point.
          RSeQC is excluded until a gene BED file is provided via --gene_bed.
========================================================================================
*/

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT  } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS  } from '../../modules/nf-core/samtools/idxstats/main'

workflow BAM_QC {

    take:
    ch_bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
    ch_gene_bed  // channel: not used yet — reserved for future RSeQC integration

    main:

    ch_multiqc_files = Channel.empty()

    // ── samtools flagstat — alignment summary stats ───────────────────────────
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_multiqc_files = ch_multiqc_files.mix(
        SAMTOOLS_FLAGSTAT.out.flagstat.collect{ it[1] }.ifEmpty([])
    )

    // ── samtools idxstats — per-chromosome read counts ───────────────────────
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_multiqc_files = ch_multiqc_files.mix(
        SAMTOOLS_IDXSTATS.out.idxstats.collect{ it[1] }.ifEmpty([])
    )

    emit:
    multiqc_files = ch_multiqc_files
}
