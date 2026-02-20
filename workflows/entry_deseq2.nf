/*
========================================================================================
    ENTRY POINT: DESEQ2
    Re-runs visualizations and the HTML report from a saved DESeq2 RDS object.
    Useful for: adjusting thresholds, adding genes of interest, tweaking plot aesthetics.

    USAGE:
        nextflow run main.nf -profile docker \\
            --entry_point deseq2 \\
            --deseq2_rds ./results/deseq2/deseq2_object.rds \\
            --input samplesheet.csv \\
            --outdir ./results_replot \\
            --lfc_threshold 0.5 \\
            --genes_of_interest "IGF2BP3,MYC,BCL2"
========================================================================================
*/

nextflow.enable.dsl = 2

include { DESEQ2_REPLOT         } from '../modules/local/deseq2_replot/main'
include { GENERATE_REPORT       } from '../modules/local/generate_report/main'

workflow PIPELINE_FROM_DESEQ2 {

    ch_versions = Channel.empty()

    // Wrap the RDS file in a channel
    ch_rds = Channel.value(file(params.deseq2_rds))

    //
    // Re-extract results + regenerate all plots from saved RDS
    //
    DESEQ2_REPLOT(
        ch_rds,
        file(params.input),
        params.ko_condition,
        params.control_condition,
        params.lfc_threshold,
        params.pval_threshold,
        params.top_de_genes,
        params.genes_of_interest ?: ''
    )
    ch_versions = ch_versions.mix(DESEQ2_REPLOT.out.versions)

    //
    // Re-generate HTML report
    //
    if (!params.skip_report) {
        GENERATE_REPORT(
            DESEQ2_REPLOT.out.results_all,
            DESEQ2_REPLOT.out.results_sig,
            DESEQ2_REPLOT.out.ko_stats,
            DESEQ2_REPLOT.out.ko_plot,
            DESEQ2_REPLOT.out.plots,
            file(params.input),
            ch_versions.collect()
        )
    }
}
