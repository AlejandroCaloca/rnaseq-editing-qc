/*
========================================================================================
    LOCAL MODULE: DESEQ2_REPLOT
    Re-extract results and regenerate all plots from a saved DESeq2 RDS object.
    Used by --entry_point deseq2 to adjust thresholds or aesthetics without re-running DE.
========================================================================================
*/

process DESEQ2_REPLOT {
    label 'process_low'
    publishDir "${params.outdir}/deseq2", mode: params.publish_dir_mode

    conda "conda-forge::r-base=4.3.0 bioconda::bioconductor-deseq2=1.40.0 conda-forge::r-ggplot2=3.4.2 conda-forge::r-pheatmap conda-forge::r-ggrepel conda-forge::r-openxlsx"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.40.0--r43hf17093f_0' :
        'biocontainers/bioconductor-deseq2:1.40.0--r43hf17093f_0' }"

    input:
    path rds_file
    path samplesheet
    val  ko_condition
    val  control_condition
    val  lfc_threshold
    val  pval_threshold
    val  top_n_genes
    val  genes_of_interest

    output:
    path "deseq2_results_all.csv",  emit: results_all
    path "deseq2_results_sig.csv",  emit: results_sig
    path "ko_efficiency_stats.csv", emit: ko_stats
    path "ko_efficiency_plot.png",  emit: ko_plot
    path "plots/",                  emit: plots
    path "versions.yml",            emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(openxlsx)

    dir.create("plots", showWarnings = FALSE)

    # Load saved DESeq2 object
    cat("Loading DESeq2 object from: ${rds_file}\\n")
    dds <- readRDS("${rds_file}")

    # Re-extract results with new thresholds
    res <- results(dds,
                   contrast = c("condition", "${ko_condition}", "${control_condition}"),
                   alpha    = ${pval_threshold})

    res_df  <- as.data.frame(res)
    res_df\$gene_id <- rownames(res_df)
    res_df  <- res_df[order(res_df\$padj, na.last = TRUE), ]
    res_sig <- res_df[!is.na(res_df\$padj) &
                       res_df\$padj < ${pval_threshold} &
                       abs(res_df\$log2FoldChange) >= ${lfc_threshold}, ]

    write.csv(res_df,  "deseq2_results_all.csv", row.names = FALSE)
    write.csv(res_sig, "deseq2_results_sig.csv", row.names = FALSE)
    write.xlsx(res_df,  "deseq2_results_all.xlsx", rowNames = FALSE)
    write.xlsx(res_sig, "deseq2_results_sig.xlsx", rowNames = FALSE)

    cat("Significant genes:", nrow(res_sig), "\\n")
    cat("  Up:  ", sum(res_sig\$log2FoldChange > 0), "\\n")
    cat("  Down:", sum(res_sig\$log2FoldChange < 0), "\\n")

    # Variance-stabilised counts for plots
    vsd <- vst(dds, blind = FALSE)

    # PCA
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
    pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text_repel(size = 3) +
        scale_color_manual(values = c("${ko_condition}" = "#E74C3C",
                                       "${control_condition}" = "#3498DB")) +
        labs(x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)"),
             title = "PCA: Sample Clustering", color = "Condition") +
        theme_bw()
    ggsave("plots/pca_plot.pdf", pca_plot, width = 8, height = 6)
    ggsave("plots/pca_plot.png", pca_plot, width = 8, height = 6, dpi = 150)

    # Volcano
    pd <- res_df[!is.na(res_df\$padj), ]
    pd\$sig <- "NS"
    pd\$sig[pd\$padj < ${pval_threshold} & pd\$log2FoldChange >= ${lfc_threshold}]  <- "Up"
    pd\$sig[pd\$padj < ${pval_threshold} & pd\$log2FoldChange <= -${lfc_threshold}] <- "Down"
    pd\$label <- ifelse(pd\$sig != "NS" & rank(pd\$padj) <= 20, pd\$gene_id, "")

    volcano <- ggplot(pd, aes(log2FoldChange, -log10(padj), color = sig, label = label)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(size = 2.5, max.overlaps = 15) +
        geom_vline(xintercept = c(-${lfc_threshold}, ${lfc_threshold}), linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(${pval_threshold}), linetype = "dashed", color = "grey50") +
        scale_color_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "grey70")) +
        labs(title = "Volcano Plot: ${ko_condition} vs ${control_condition}",
             x = "log2 Fold Change", y = "-log10(padj)", color = "") +
        theme_bw()
    ggsave("plots/volcano_plot.pdf", volcano, width = 10, height = 8)
    ggsave("plots/volcano_plot.png", volcano, width = 10, height = 8, dpi = 150)

    # MA plot
    ma_d <- res_df[!is.na(res_df\$padj), ]
    ma_d\$significant <- ma_d\$padj < ${pval_threshold}
    ma <- ggplot(ma_d, aes(log10(baseMean + 1), log2FoldChange, color = significant)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = c(-${lfc_threshold}, ${lfc_threshold}), linetype = "dashed", color = "grey40") +
        scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "grey70"),
                           labels = c("TRUE" = "Significant", "FALSE" = "Not Sig.")) +
        labs(title = "MA Plot: ${ko_condition} vs ${control_condition}",
             x = "log10(Mean Expression + 1)", y = "log2 Fold Change") +
        theme_bw()
    ggsave("plots/ma_plot.pdf", ma, width = 10, height = 6)
    ggsave("plots/ma_plot.png", ma, width = 10, height = 6, dpi = 150)

    # Heatmap
    top_genes <- head(res_sig\$gene_id, ${top_n_genes})
    if (length(top_genes) >= 2) {
        mat <- assay(vsd)[top_genes, ]
        mat <- mat - rowMeans(mat)
        annot_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
        pdf("plots/heatmap_top_genes.pdf", width = 10, height = 14)
        pheatmap(mat, annotation_col = annot_col,
                 main = paste("Top", length(top_genes), "DE Genes"),
                 color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100))
        dev.off()
        png("plots/heatmap_top_genes.png", width = 1000, height = 1400, res = 100)
        pheatmap(mat, annotation_col = annot_col,
                 main = paste("Top", length(top_genes), "DE Genes"),
                 color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100))
        dev.off()
    }

    # Genes of interest barplots
    goi_str <- "${genes_of_interest}"
    if (nchar(goi_str) > 0) {
        goi <- trimws(strsplit(goi_str, ",")[[1]])
        goi_present <- goi[goi %in% rownames(vsd)]
        if (length(goi_present) > 0) {
            goi_mat <- as.data.frame(t(assay(vsd)[goi_present, , drop = FALSE]))
            goi_mat\$sample <- rownames(goi_mat)
            goi_mat\$condition <- colData(vsd)\$condition
            for (gene in goi_present) {
                p <- ggplot(goi_mat, aes(x = sample, y = .data[[gene]], fill = condition)) +
                    geom_bar(stat = "identity") +
                    scale_fill_manual(values = c("${ko_condition}" = "#E74C3C",
                                                  "${control_condition}" = "#3498DB")) +
                    labs(title = paste(gene, "- Normalised Expression"),
                         y = "VST counts", x = "") +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(paste0("plots/", gene, "_expression.png"), p, width = 7, height = 5, dpi = 150)
            }
        }
    }

    # Stub KO stats (re-extracted from dds object if available)
    ko_gene_counts <- tryCatch({
        counts(dds, normalized = TRUE)["${params.ko_gene}", ]
    }, error = function(e) NULL)

    if (!is.null(ko_gene_counts)) {
        ko_mean  <- mean(ko_gene_counts[colData(dds)\$condition == "${ko_condition}"])
        nt_mean  <- mean(ko_gene_counts[colData(dds)\$condition == "${control_condition}"])
        efficiency <- 1 - (ko_mean / nt_mean)
        stats_df <- data.frame(ko_gene = "${params.ko_gene}", ko_mean_count = ko_mean,
                                control_mean_count = nt_mean, ko_efficiency = efficiency,
                                threshold = ${params.ko_efficiency_threshold},
                                passed_qc = efficiency >= ${params.ko_efficiency_threshold})
        write.csv(stats_df, "ko_efficiency_stats.csv", row.names = FALSE)
        # Barplot
        plot_df <- data.frame(
            condition = c("${ko_condition}", "${control_condition}"),
            mean_count = c(ko_mean, nt_mean)
        )
        p_ko <- ggplot(plot_df, aes(condition, mean_count, fill = condition)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("${ko_condition}" = "#E74C3C",
                                          "${control_condition}" = "#3498DB")) +
            labs(title = "${params.ko_gene} Mean Normalised Expression",
                 y = "Normalised counts", x = "") + theme_bw()
        ggsave("ko_efficiency_plot.png", p_ko, width = 6, height = 5, dpi = 150)
    } else {
        write.csv(data.frame(note = "${params.ko_gene} not found in dds rownames"),
                  "ko_efficiency_stats.csv", row.names = FALSE)
        png("ko_efficiency_plot.png", width = 400, height = 300)
        plot.new(); text(0.5, 0.5, "${params.ko_gene} not found", cex = 1.5)
        dev.off()
    }

    writeLines(
        c('versions:', paste0('    "', Sys.getenv("NXF_TASK_PROCESS"), '":'),
          paste0('        r-base: "', R.version\$major, ".", R.version\$minor, '"'),
          paste0('        bioconductor-deseq2: "', packageVersion("DESeq2"), '"')),
        "versions.yml"
    )
    """
}
