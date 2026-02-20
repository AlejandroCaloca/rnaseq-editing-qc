/*
========================================================================================
    SUBWORKFLOW: BAM_QC
    Quality control of pre-aligned BAM files
    Runs: samtools flagstat, samtools idxstats, samtools stats, RSeQC (optional)
========================================================================================
*/

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT     } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS     } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_STATS        } from '../../modules/nf-core/samtools/stats/main'
include { RSEQC_INFEREXPERIMENT } from '../../modules/nf-core/rseqc/inferexperiment/main'
include { RSEQC_READDISTRIBUTION} from '../../modules/nf-core/rseqc/readdistribution/main'

workflow BAM_QC {

    take:
    ch_bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
    ch_gene_bed  // path: gene BED file (for RSeQC)

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // samtools flagstat — alignment summary stats
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{ it[1] }.ifEmpty([]))
    ch_versions      = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    // samtools idxstats — per-chromosome read counts
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect{ it[1] }.ifEmpty([]))
    ch_versions      = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    // samtools stats — comprehensive stats
    SAMTOOLS_STATS(ch_bam_bai, [])
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{ it[1] }.ifEmpty([]))

    // RSeQC: infer strandedness — important for confirming HTSeq --stranded setting
    if (ch_gene_bed && !params.skip_rseqc) {
        RSEQC_INFEREXPERIMENT(ch_bam_bai, ch_gene_bed)
        ch_multiqc_files = ch_multiqc_files.mix(
            RSEQC_INFEREXPERIMENT.out.txt.collect{ it[1] }.ifEmpty([])
        )
        ch_versions = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())

        // RSeQC: read distribution across genomic features
        RSEQC_READDISTRIBUTION(ch_bam_bai, ch_gene_bed)
        ch_multiqc_files = ch_multiqc_files.mix(
            RSEQC_READDISTRIBUTION.out.txt.collect{ it[1] }.ifEmpty([])
        )
    }

    emit:
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
