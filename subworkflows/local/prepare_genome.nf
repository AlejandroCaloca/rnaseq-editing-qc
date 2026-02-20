/*
========================================================================================
    SUBWORKFLOW: PREPARE_GENOME
    Prepares genome indices for alignment and quantification
========================================================================================
*/

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE   } from '../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD          } from '../../modules/nf-core/hisat2/build/main'
include { SALMON_INDEX          } from '../../modules/nf-core/salmon/index/main'

workflow PREPARE_GENOME {

    main:
    ch_versions = Channel.empty()

    // Set genome files from params or iGenomes
    ch_fasta            = params.fasta            ? file(params.fasta)            : Channel.empty()
    ch_gtf              = params.gtf              ? file(params.gtf)              : Channel.empty()
    ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta) : Channel.empty()

    // Build STAR index if not provided
    ch_star_index = Channel.empty()
    if (params.aligner == 'star') {
        if (params.star_index) {
            ch_star_index = Channel.value(file(params.star_index))
        } else {
            STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
            ch_star_index = STAR_GENOMEGENERATE.out.index
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    // Build HISAT2 index if not provided
    ch_hisat2_index  = Channel.empty()
    ch_splicesites   = Channel.empty()
    if (params.aligner == 'hisat2') {
        if (params.hisat2_index) {
            ch_hisat2_index = Channel.value(file(params.hisat2_index))
        } else {
            HISAT2_EXTRACTSPLICESITES(ch_gtf)
            ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt
            HISAT2_BUILD(ch_fasta, ch_gtf, ch_splicesites)
            ch_hisat2_index = HISAT2_BUILD.out.index
            ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    // Build Salmon index
    ch_salmon_index = Channel.empty()
    if (params.quantification_method == 'salmon' || params.quantification_method == 'both') {
        if (params.salmon_index) {
            ch_salmon_index = Channel.value(file(params.salmon_index))
        } else {
            SALMON_INDEX(ch_fasta, ch_transcript_fasta)
            ch_salmon_index = SALMON_INDEX.out.index
            ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    emit:
    fasta            = ch_fasta
    gtf              = ch_gtf
    transcript_fasta = ch_transcript_fasta
    star_index       = ch_star_index
    hisat2_index     = ch_hisat2_index
    splicesites      = ch_splicesites
    salmon_index     = ch_salmon_index
    versions         = ch_versions
}
