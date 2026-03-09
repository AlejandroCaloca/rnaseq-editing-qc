/*
========================================================================================
    LOCAL MODULE: DESEQ2_ANALYSIS
    Differential expression analysis for KO vs Control
========================================================================================
*/

process DESEQ2_ANALYSIS {
    label 'process_medium'
    publishDir "${params.outdir}/deseq2", mode: params.publish_dir_mode

   conda "conda-forge::r-base=4.3.1 bioconda::bioconductor-deseq2=1.40.2 conda-forge::r-ggplot2 conda-forge::r-pheatmap conda-forge::r-ggrepel conda-forge::r-openxlsx conda-forge::r-dplyr"
        container null

  //  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
  //   'docker://rocker/tidyverse:4.3.1' : 'rocker/tidyverse:4.3.1' }"  

   // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   //     'docker://bioconductor/bioconductor_docker:RELEASE_3_18' :
   //     'bioconductor/bioconductor_docker:RELEASE_3_18' }"
    
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.42.0--r43hdfd78af_0' :
    //    'quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hdfd78af_0' }"



    input:
    path count_files
    path samplesheet
    val ko_condition
    val control_condition
    val lfc_threshold
    val pval_threshold

    output:
    path "deseq2_results_all.csv",      emit: results_all
    path "deseq2_results_all.xlsx",     emit: results_all_xlsx
    path "deseq2_results_sig.csv",      emit: results_sig
    path "deseq2_results_sig.xlsx",     emit: results_sig_xlsx
    path "deseq2_object.rds",           emit: rds
    path "plots/",                      emit: plots
    path "deseq2_summary_mqc.json",     emit: multiqc_data
    path "versions.yml",                emit: versions

    script:
    """
    Rscript - <<'EOF'

    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(openxlsx)
    library(dplyr)

    dir.create("plots", showWarnings = FALSE)

    # ---------------------------------------------------------------------------------
    # 1. Build count matrix from HTSeq files
    # ---------------------------------------------------------------------------------
    count_files <- list.files(".", pattern = "\\\\.counts\\\\.txt\$", full.names = TRUE)
    
    read_htseq <- function(f) {
        df <- read.table(f, header = FALSE, sep = "\\t",
                         col.names = c("gene_id", "gene_name", "count"))
        df <- df[!grepl("^__", df\$gene_id), ]
        rownames(df) <- df\$gene_id
        setNames(df["count"], gsub("\\\\.counts\\\\.txt", "", basename(f)))
    }
    
    count_list <- lapply(count_files, read_htseq)
    count_matrix <- do.call(cbind, count_list)
    count_matrix <- as.matrix(count_matrix)
    mode(count_matrix) <- "integer"
    
    # ---------------------------------------------------------------------------------
    # 2. Build colData from samplesheet
    # ---------------------------------------------------------------------------------
    samplesheet <- read.csv("${samplesheet}")
    coldata <- samplesheet[match(colnames(count_matrix), samplesheet\$sample), ]
    rownames(coldata) <- coldata\$sample
    coldata\$condition <- factor(coldata\$condition, 
                                  levels = c("${control_condition}", "${ko_condition}"))
    
    # ---------------------------------------------------------------------------------
    # 3. Run DESeq2
    # ---------------------------------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData   = coldata,
        design    = ~ condition
    )
    
    # Filter low-count genes
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds  <- dds[keep, ]
    
    dds <- DESeq(dds)
    saveRDS(dds, "deseq2_object.rds")
    
    # ---------------------------------------------------------------------------------
    # 4. Extract results
    # ---------------------------------------------------------------------------------
    res <- results(dds,
                   contrast  = c("condition", "${ko_condition}", "${control_condition}"),
                   alpha     = ${pval_threshold},
                   lfcThreshold = 0)  # Report all; filter below
    
    res_df <- as.data.frame(res)
    res_df\$gene_id <- rownames(res_df)
    res_df <- res_df[order(res_df\$padj, na.last = TRUE), ]
    
    # Save all results
    write.csv(res_df, "deseq2_results_all.csv", row.names = FALSE)
    write.xlsx(res_df, "deseq2_results_all.xlsx", rowNames = FALSE)
    
    # Significant results
    res_sig <- res_df[!is.na(res_df\$padj) &
                       res_df\$padj < ${pval_threshold} &
                       abs(res_df\$log2FoldChange) >= ${lfc_threshold}, ]
    
    write.csv(res_sig, "deseq2_results_sig.csv", row.names = FALSE)
    write.xlsx(res_sig, "deseq2_results_sig.xlsx", rowNames = FALSE)
    
    n_up   <- sum(res_sig\$log2FoldChange > 0)
    n_down <- sum(res_sig\$log2FoldChange < 0)
    
    cat("\\n========================================\\n")
    cat("DESeq2 Analysis Summary\\n")
    cat("========================================\\n")
    cat("Contrast: ${ko_condition} vs ${control_condition}\\n")
    cat("Total genes tested: ", nrow(res_df), "\\n")
    cat("Significant (padj<${pval_threshold}, |LFC|>${lfc_threshold}): ", nrow(res_sig), "\\n")
    cat("  Upregulated:  ", n_up, "\\n")
    cat("  Downregulated:", n_down, "\\n")
    cat("========================================\\n\\n")
    
    # ---------------------------------------------------------------------------------
    # 5. PCA Plot
    # ---------------------------------------------------------------------------------
    vsd <- vst(dds, blind = FALSE)
    
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
    
    pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4, alpha = 0.9) +
        geom_text_repel(size = 3, max.overlaps = 20) +
        scale_color_manual(values = c("${ko_condition}" = "#E74C3C",
                                       "${control_condition}" = "#3498DB")) +
        labs(
            title = "PCA: Sample Clustering",
            x = paste0("PC1 (", pct_var[1], "%)"),
            y = paste0("PC2 (", pct_var[2], "%)"),
            color = "Condition"
        ) +
        theme_bw()
    
    ggsave("plots/pca_plot.pdf", pca_plot, width = 8, height = 6)
    ggsave("plots/pca_plot.png", pca_plot, width = 8, height = 6, dpi = 150)
    
    # ---------------------------------------------------------------------------------
    # 6. Volcano Plot
    # ---------------------------------------------------------------------------------
    plot_data <- res_df[!is.na(res_df\$padj), ]
    plot_data\$significance <- "Not Significant"
    plot_data\$significance[plot_data\$padj < ${pval_threshold} & 
                             plot_data\$log2FoldChange >= ${lfc_threshold}]  <- "Upregulated"
    plot_data\$significance[plot_data\$padj < ${pval_threshold} & 
                             plot_data\$log2FoldChange <= -${lfc_threshold}] <- "Downregulated"
    
    plot_data\$label <- ifelse(
        plot_data\$significance != "Not Significant" & 
        rank(plot_data\$padj) <= 20,
        plot_data\$gene_id, ""
    )
    
    volcano <- ggplot(plot_data, aes(log2FoldChange, -log10(padj), 
                                     color = significance, label = label)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(size = 2.5, max.overlaps = 15) +
        geom_vline(xintercept = c(-${lfc_threshold}, ${lfc_threshold}), 
                   linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(${pval_threshold}), 
                   linetype = "dashed", color = "grey50") +
        scale_color_manual(values = c("Upregulated" = "#E74C3C",
                                       "Downregulated" = "#3498DB",
                                       "Not Significant" = "grey70")) +
        labs(
            title = "Volcano Plot: ${ko_condition} vs ${control_condition}",
            x = "log2 Fold Change",
            y = "-log10(adjusted p-value)",
            color = "Significance"
        ) +
        theme_bw()
    
    ggsave("plots/volcano_plot.pdf", volcano, width = 10, height = 8)
    ggsave("plots/volcano_plot.png", volcano, width = 10, height = 8, dpi = 150)
    
    # ---------------------------------------------------------------------------------
    # 7. MA Plot
    # ---------------------------------------------------------------------------------
    ma_data <- res_df[!is.na(res_df\$padj), ]
    ma_data\$significant <- ma_data\$padj < ${pval_threshold}
    
    ma_plot <- ggplot(ma_data, aes(log10(baseMean + 1), log2FoldChange, 
                                    color = significant)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_hline(yintercept = 0, color = "black") +
        geom_hline(yintercept = c(-${lfc_threshold}, ${lfc_threshold}), 
                   linetype = "dashed", color = "grey40") +
        scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "grey70"),
                           labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
        labs(
            title = "MA Plot: ${ko_condition} vs ${control_condition}",
            x = "log10(Mean Expression + 1)",
            y = "log2 Fold Change",
            color = ""
        ) +
        theme_bw()
    
    ggsave("plots/ma_plot.pdf", ma_plot, width = 10, height = 6)
    ggsave("plots/ma_plot.png", ma_plot, width = 10, height = 6, dpi = 150)
    
    # ---------------------------------------------------------------------------------
    # 8. Heatmap (top DE genes)
    # ---------------------------------------------------------------------------------
    top_genes <- head(res_sig\$gene_id, ${params.top_de_genes})
    
    if (length(top_genes) >= 2) {
        mat <- assay(vsd)[top_genes, ]
        mat <- mat - rowMeans(mat)
        
        annot_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
        annot_colors <- list(condition = c("${ko_condition}" = "#E74C3C",
                                            "${control_condition}" = "#3498DB"))
        
        pdf("plots/heatmap_top_genes.pdf", width = 10, height = 14)
        pheatmap(mat,
                 annotation_col  = annot_col,
                 annotation_colors = annot_colors,
                 cluster_rows    = ${params.heatmap_cluster_rows ? "TRUE" : "FALSE"},
                 cluster_cols    = ${params.heatmap_cluster_cols ? "TRUE" : "FALSE"},
                 show_rownames   = length(top_genes) <= 100,
                 main            = paste0("Top ", length(top_genes), " DE Genes"),
                 color           = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100))
        dev.off()
        
        png("plots/heatmap_top_genes.png", width = 1000, height = 1400, res = 100)
        pheatmap(mat,
                 annotation_col  = annot_col,
                 annotation_colors = annot_colors,
                 cluster_rows    = ${params.heatmap_cluster_rows ? "TRUE" : "FALSE"},
                 cluster_cols    = ${params.heatmap_cluster_cols ? "TRUE" : "FALSE"},
                 show_rownames   = length(top_genes) <= 100,
                 main            = paste0("Top ", length(top_genes), " DE Genes"),
                 color           = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100))
        dev.off()
    }
    
    # ---------------------------------------------------------------------------------
    # 9. MultiQC JSON summary
    # ---------------------------------------------------------------------------------
    mqc_json <- paste0('{
        "id": "deseq2_summary",
        "section_name": "DESeq2 Differential Expression",
        "description": "Summary of differential expression analysis (${ko_condition} vs ${control_condition})",
        "plot_type": "generalstats",
        "pconfig": [
            {"Total Tested": {"title": "Total Tested", "format": "{:,.0f}"}},
            {"Significant": {"title": "Significant", "format": "{:,.0f}"}},
            {"Upregulated": {"title": "Upregulated", "format": "{:,.0f}"}},
            {"Downregulated": {"title": "Downregulated", "format": "{:,.0f}"}}
        ],
        "data": {
            "${ko_condition}_vs_${control_condition}": {
                "Total Tested": ', nrow(res_df), ',
                "Significant": ', nrow(res_sig), ',
                "Upregulated": ', n_up, ',
                "Downregulated": ', n_down, '
            }
        }
    }')
    writeLines(mqc_json, "deseq2_summary_mqc.json")

    writeLines(
        c('versions:', 
          paste0('    "', Sys.getenv("NXF_TASK_PROCESS"), '":'),
          paste0('        r-base: "', R.version\$major, ".", R.version\$minor, '"'),
          paste0('        bioconductor-deseq2: "', packageVersion("DESeq2"), '"')),
        "versions.yml"
    )

    EOF
    
    """
}
