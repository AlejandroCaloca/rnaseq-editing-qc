/*
========================================================================================
    LOCAL MODULE: GENERATE_REPORT
    Creates the final HTML report via RMarkdown
========================================================================================
*/

process GENERATE_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    conda "conda-forge::r-base=4.3.0 conda-forge::r-rmarkdown conda-forge::r-knitr conda-forge::r-kableextra conda-forge::r-dt conda-forge::r-ggplot2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.40.0--r43hf17093f_0' :
        'quay.io/biocontainers/bioconductor-deseq2:1.40.0--r43hf17093f_0' }"

    input:
    path results_all
    path results_sig
    path ko_stats
    path ko_plot
    path plots_dir
    path samplesheet
    path versions

    output:
    path "*.html",  emit: report
    path "*.pdf",   emit: report_pdf, optional: true

    script:
    def report_title = params.report_title ?: 'RNA-seq Analysis Report'
    """
    # Copy report template and render
    cp ${workflow.projectDir}/assets/report_template.Rmd report.Rmd

    Rscript -e "
    rmarkdown::render(
        'report.Rmd',
        output_file  = '${params.ko_condition}_vs_${params.control_condition}_analysis_report.html',
        params = list(
            title              = '${report_title}',
            results_all        = '${results_all}',
            results_sig        = '${results_sig}',
            ko_stats           = '${ko_stats}',
            ko_condition       = '${params.ko_condition}',
            control_condition  = '${params.control_condition}',
            ko_gene            = '${params.ko_gene}',
            lfc_threshold      = ${params.lfc_threshold},
            pval_threshold     = ${params.pval_threshold},
            samplesheet        = '${samplesheet}'
        )
    )
    "
    """
}
