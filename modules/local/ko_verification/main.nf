/*
========================================================================================
    LOCAL MODULE: KO_VERIFICATION
    Checks IGF2BP3 knockout efficiency across all samples
========================================================================================
*/

process KO_VERIFICATION {
    label 'process_low'
    publishDir "${params.outdir}/ko_verification", mode: params.publish_dir_mode

    conda "conda-forge::r-base=4.3.0 bioconda::bioconductor-deseq2=1.40.0 conda-forge::r-ggplot2=3.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.40.0--r43hf17093f_0' :
        'quay.io/biocontainers/bioconductor-deseq2:1.40.0--r43hf17093f_0' }"

    input:
    path count_files
    val ko_gene
    val ko_condition
    val control_condition
    val threshold

    output:
    path "ko_efficiency_stats.csv",     emit: ko_stats
    path "ko_efficiency_plot.pdf",      emit: ko_plot
    path "ko_efficiency_plot.png",      emit: ko_plot_png
    path "ko_verification_mqc.json",    emit: multiqc_data
    path "versions.yml",                emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)

    # Read all count files
    count_files <- list.files(".", pattern = "\\\\.counts\\\\.txt\$", full.names = TRUE)
    
    # Parse counts into a matrix
    read_counts <- function(f) {
        df <- read.table(f, header = FALSE, sep = "\\t", comment.char = "#",
                         col.names = c("gene_id", "gene_name", "count"))
        df <- df[!grepl("^__", df\$gene_id), ]
        sample_name <- gsub("\\\\.counts\\\\.txt", "", basename(f))
        df\$sample <- sample_name
        return(df)
    }
    
    all_counts <- do.call(rbind, lapply(count_files, read_counts))
    
    # Extract KO gene counts
    ko_data <- all_counts[all_counts\$gene_name == "${ko_gene}" | 
                           all_counts\$gene_id == "${ko_gene}", ]
    
    if (nrow(ko_data) == 0) {
        stop("ERROR: Gene '${ko_gene}' not found in count files. Check gene name spelling.")
    }
    
    # Assign condition based on sample name
    ko_data\$condition <- ifelse(
        grepl("${ko_condition}", ko_data\$sample), "${ko_condition}", "${control_condition}"
    )
    
    # Calculate mean by condition
    mean_counts <- ko_data %>%
        group_by(condition) %>%
        summarise(mean_count = mean(count), sd_count = sd(count), .groups = 'drop')
    
    ko_mean  <- mean_counts\$mean_count[mean_counts\$condition == "${ko_condition}"]
    nt_mean  <- mean_counts\$mean_count[mean_counts\$condition == "${control_condition}"]
    
    ko_efficiency <- 1 - (ko_mean / nt_mean)
    passed        <- ko_efficiency >= ${threshold}
    
    # Write stats
    stats_df <- data.frame(
        ko_gene           = "${ko_gene}",
        ko_mean_count     = ko_mean,
        control_mean_count = nt_mean,
        ko_efficiency     = ko_efficiency,
        threshold         = ${threshold},
        passed_qc         = passed
    )
    write.csv(stats_df, "ko_efficiency_stats.csv", row.names = FALSE)
    
    # Generate plot
    p <- ggplot(ko_data, aes(x = sample, y = count, fill = condition)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("${ko_condition}" = "#E74C3C", "${control_condition}" = "#3498DB")) +
        labs(
            title    = paste0("${ko_gene} Knockout Verification"),
            subtitle = paste0("KO efficiency: ", round(ko_efficiency * 100, 1), "% | ",
                              ifelse(passed, "PASSED QC ✓", "FAILED QC ✗")),
            x        = "Sample",
            y        = "Raw Read Count",
            fill     = "Condition"
        ) +
        theme_bw() +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1),
            plot.subtitle = element_text(
                color = ifelse(passed, "darkgreen", "red"), face = "bold"
            )
        )
    
    ggsave("ko_efficiency_plot.pdf", p, width = 10, height = 6)
    ggsave("ko_efficiency_plot.png", p, width = 10, height = 6, dpi = 150)
    
    # Write MultiQC-compatible JSON
    mqc_json <- paste0('{
        "id": "ko_verification",
        "section_name": "Knockout Verification",
        "description": "IGF2BP3 expression levels confirming knockout efficiency",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "ko_bargraph",
            "title": "${ko_gene} Expression by Sample",
            "ylab": "Read Counts"
        },
        "data": {', 
        paste(sapply(1:nrow(ko_data), function(i) {
            paste0('"', ko_data\$sample[i], '": {"${ko_gene}": ', ko_data\$count[i], '}')
        }), collapse = ", "),
        '}}'
    )
    writeLines(mqc_json, "ko_verification_mqc.json")
    
    # Print status
    cat("\\n========================================\\n")
    cat("Knockout Verification Summary\\n")
    cat("========================================\\n")
    cat("Gene:            ${ko_gene}\\n")
    cat("KO mean count:  ", ko_mean, "\\n")
    cat("NT mean count:  ", nt_mean, "\\n")
    cat("KO efficiency:  ", round(ko_efficiency * 100, 1), "%\\n")
    cat("Threshold:       ${threshold}%\\n")
    cat("QC status:      ", ifelse(passed, "PASSED", "FAILED"), "\\n")
    cat("========================================\\n\\n")
    
    if (!passed) {
        warning(paste0(
            "WARNING: Knockout efficiency (", round(ko_efficiency * 100, 1), "%) ",
            "is below the ", ${threshold} * 100, "% threshold. ",
            "Check samples for incomplete editing."
        ))
    }

    writeLines(
        c('versions:', paste0('    "', Sys.getenv("NXF_TASK_PROCESS"), '":'),
          paste0('        r-base: "', R.version\$major, ".", R.version\$minor, '"')),
        "versions.yml"
    )
    """
}
