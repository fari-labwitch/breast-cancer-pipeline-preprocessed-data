# src/visualizer.R

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

#' Save a ggplot or patchwork plot robustly
#'
#' @param plot_obj The ggplot/patchwork object
#' @param filename Character path to save the plot
#' @param width Numeric width in inches
#' @param height Numeric height in inches
#' @param dpi Numeric resolution in dpi
save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 150) {
  tryCatch({
    message(sprintf("Saving plot to %s ...", filename))
    ggsave(filename = filename, plot = plot_obj, width = width, height = height, dpi = dpi)
    message("Saved successfully.")
  }, error = function(e) {
    warning(sprintf("Failed to save plot %s: %s", filename, e$message))
  })
}

#' Sanitize gost results by converting list columns (like 'parents') into character strings.
#'
#' @description
#' The `gprofiler2::gost()` result contains some columns (such as `parents`) that are stored
#' as lists of character vectors in R. For example, a single row might contain:
#'   `parents = c("GO:0008150", "GO:0050789")`
#' Standard R serialization functions like `write.csv()` or `write.table()` cannot write list
#' columns directly, causing:
#'   `Error in utils::write.table(...) : unimplemented type 'list' in 'EncodeElement'`
#' This function scans the input dataframe for columns of type list and collapses them
#' into a single comma-separated string (e.g. `"GO:0008150, GO:0050789"`), which can then
#' be safely written to CSV.
#'
#' @param df The raw gost results dataframe.
#' @return A sanitized dataframe with list columns converted to character vectors.
sanitize_gost_result <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  # Find columns that are lists
  list_cols <- sapply(df, is.list)
  
  # Convert each list column to a comma-separated string
  for (col_name in names(df)[list_cols]) {
    df[[col_name]] <- sapply(df[[col_name]], function(x) {
      if (is.null(x)) {
        return("")
      } else {
        return(paste(x, collapse = ", "))
      }
    })
  }
  return(df)
}

#' Generate and Save Assignment Plots
#'
#' @param se A processed Seurat object
#' @param results_dir Character path to the results directory
generate_assignment_plots <- function(se, results_dir = "results") {
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Ensure we have active devices configured for headless environment
  # Set pdf(NULL) to prevent opening GUI windows
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  
  # 1. Initial UMAP Plot
  message("Generating initial UMAP plot...")
  # Keep track of current default identity
  orig_ident <- Idents(se)
  
  p_init <- DimPlot(se, reduction = "umap") + 
    labs(title = "UMAP Dimensional Projection (Initial)") +
    theme_minimal()
  save_plot(p_init, file.path(results_dir, "umap_initial.png"), width = 7, height = 6)
  
  # 2. UMAP of clusters at resolution 0.6
  message("Generating UMAP of clusters (res 0.6)...")
  p_res06 <- DimPlot(se, reduction = "umap", group.by = "SCT_snn_res.0.6") +
    labs(title = "Clustering at Resolution 0.6") +
    theme_minimal()
  save_plot(p_res06, file.path(results_dir, "umap_clusters_res_0.6.png"), width = 7, height = 6)
  
  # 3. UMAPs of clusters at resolution 0.3 and 1.2
  message("Generating UMAP of clusters (res 0.3 & 1.2)...")
  p_res03 <- DimPlot(se, reduction = "umap", group.by = "SCT_snn_res.0.3") +
    labs(title = "Clustering at Resolution 0.3") +
    theme_minimal()
  p_res12 <- DimPlot(se, reduction = "umap", group.by = "SCT_snn_res.1.2") +
    labs(title = "Clustering at Resolution 1.2") +
    theme_minimal()
  
  save_plot(p_res03, file.path(results_dir, "umap_clusters_res_0.3.png"), width = 7, height = 6)
  save_plot(p_res12, file.path(results_dir, "umap_clusters_res_1.2.png"), width = 7, height = 6)
  
  p_comb_res <- p_res03 + p_res12
  save_plot(p_comb_res, file.path(results_dir, "umap_clusters_res_0.3_1.2.png"), width = 12, height = 6)
  
  # 4. UMAP showing CellType annotations
  message("Generating UMAP of cell types...")
  p_celltypes <- DimPlot(se, group.by = "CellType", label = TRUE, repel = TRUE) +
    labs(title = "Cell Type Annotations") +
    theme_minimal() +
    theme(legend.position = "right")
  save_plot(p_celltypes, file.path(results_dir, "umap_celltypes.png"), width = 9, height = 7)
  
  # 5. Marker Expression Plots (MS4A1, KIT, CD14)
  message("Generating expression plots for canonical markers (MS4A1, KIT, CD14)...")
  p_feat_markers <- FeaturePlot(se, features = c("MS4A1", "KIT", "CD14"), reduction = "umap", ncol = 3) +
    plot_annotation(title = "Marker Expression Profiles (UMAP)")
  save_plot(p_feat_markers, file.path(results_dir, "featureplot_markers_MS4A1_KIT_CD14.png"), width = 12, height = 4)
  
  p_vln_markers <- VlnPlot(se, features = c("MS4A1", "KIT", "CD14"), group.by = "CellType") +
    plot_annotation(title = "Marker Expression by Cell Type (Violin)")
  save_plot(p_vln_markers, file.path(results_dir, "vlnplot_markers_MS4A1_KIT_CD14.png"), width = 10, height = 8)
  
  # 6. Top 1 marker per cluster in FeaturePlot (based on res 0.6 clustering)
  markers <- se@misc$markers
  if (!is.null(markers) && nrow(markers) > 0) {
    message("Generating top 1 marker per cluster FeaturePlot...")
    top1_markers <- markers %>%
      group_by(cluster) %>%
      top_n(n = 1, wt = avg_log2FC) %>%
      pull(gene) %>%
      unique()
    
    # Grid size configuration
    num_cols <- 4
    p_top1 <- FeaturePlot(se, features = top1_markers, reduction = "umap", ncol = num_cols) +
      plot_annotation(title = "Top 1 Marker Gene per Cluster (UMAP)")
    save_plot(p_top1, file.path(results_dir, "featureplot_top1_markers_per_cluster.png"), 
              width = 14, height = 15)
    
    # 7. Heatmap of top 6 markers per cluster
    message("Generating top 6 markers Heatmap...")
    top6_markers <- markers %>%
      group_by(cluster) %>%
      top_n(n = 6, wt = avg_log2FC)
    
    # Set identity to SCT_snn_res.0.6 to align heatmap headers
    se <- SetIdent(se, value = "SCT_snn_res.0.6")
    p_heatmap <- DoHeatmap(se, features = top6_markers$gene, group.by = "CellType") +
      labs(title = "Top 6 Markers per Cell Type Heatmap")
    save_plot(p_heatmap, file.path(results_dir, "heatmap_top6_markers.png"), width = 12, height = 14)
  }
  
  # Restore original identity
  Idents(se) <- orig_ident
  
  # 8. Cell Cycle Violin and Feature Plots
  message("Generating Cell Cycle scoring plots...")
  p_vln_cc <- VlnPlot(se, features = c("S.Score", "G2M.Score"), group.by = "CellType") +
    plot_annotation(title = "Cell Cycle Phase Scores by Cell Type")
  save_plot(p_vln_cc, file.path(results_dir, "vlnplot_cellcycle_scores.png"), width = 12, height = 6)
  
  p_s <- FeaturePlot(se, features = "S.Score", order = TRUE, cols = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))) +
    theme_minimal()
  p_g2m <- FeaturePlot(se, features = "G2M.Score", order = TRUE, cols = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))) +
    theme_minimal()
  p_cc_feat <- p_s - p_g2m + plot_annotation(title = "Cell Cycle S.Score vs G2M.Score")
  save_plot(p_cc_feat, file.path(results_dir, "featureplot_cellcycle_scores.png"), width = 12, height = 5)
  
  # 9. Q12 Combined Patchwork Plot
  message("Generating combined patchwork plot (Q12)...")
  p1 <- FeaturePlot(se, features = "CD3E", reduction = "umap") + 
    labs(title = "T-cell Marker (CD3E) Expression") + theme_minimal()
  p2 <- DimPlot(se, reduction = "umap", group.by = "Phase") + 
    labs(title = "Cell Cycle Phase") + theme_minimal()
  p3 <- VlnPlot(se, features = "MKI67", group.by = "CellType", pt.size = 0.1) +
    labs(title = "Proliferation Marker (MKI67) Expression") + theme_minimal()
  
  combined_q12 <- (p1 | p2) / p3 + plot_annotation(title = "Proliferation and Cell Cycle State Overview")
  save_plot(combined_q12, file.path(results_dir, "patchwork_q12.png"), width = 12, height = 10)
  
  # 10. Proportions of Cell Cycle Phases by Cell Type Bar Chart
  message("Generating cell cycle phase proportion bar chart...")
  ggm <- as.data.frame(table(se@meta.data[, c("Phase", "CellType")]))
  ggm <- ggm %>%
    group_by(CellType) %>%
    mutate(proportion = Freq / sum(Freq))
  
  p_bar <- ggplot(ggm, aes(CellType, proportion, fill = Phase)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Cell Cycle Phase Proportions by Cell Type", x = "Cell Type", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1), 
          panel.grid = element_line(color = "darkgrey", linetype = "longdash")) +
    scale_fill_brewer(palette = "YlGnBu")
  save_plot(p_bar, file.path(results_dir, "barchart_cellcycle_proportions.png"), width = 10, height = 6)
  
  # 11. Functional Enrichment Analysis (FEA) results exporter
  message("Exporting FEA tables...")
  summary_file <- file.path(results_dir, "fea_summary.md")
  cat("# Functional Enrichment Analysis (FEA) Summary\n\n", file = summary_file)
  
  # Cancer/Epithelial Cell FEA (Cluster 12)
  fea_cancer <- se@misc$fea_cancer
  if (!is.null(fea_cancer) && !is.null(fea_cancer$result) && nrow(fea_cancer$result) > 0) {
    write.csv(sanitize_gost_result(fea_cancer$result), file.path(results_dir, "fea_cancer_epithelial.csv"), row.names = FALSE)
    cat("## Cancer/Epithelial Cells (Cluster 12) Enrichment\n\n", file = summary_file, append = TRUE)
    top_cancer <- head(fea_cancer$result %>% select(term_id, term_name, p_value) %>% arrange(p_value), 10)
    cat("Top 10 enriched Gene Ontology (Biological Process) terms:\n\n", file = summary_file, append = TRUE)
    cat(knitr::kable(top_cancer), sep = "\n", file = summary_file, append = TRUE)
    cat("\n\n", file = summary_file, append = TRUE)
  } else {
    cat("## Cancer/Epithelial Cells (Cluster 12) Enrichment\n\nNo enrichment results available (check log for details).\n\n", 
        file = summary_file, append = TRUE)
  }
  
  # T-cell FEA (Cluster 2)
  fea_tcell <- se@misc$fea_tcell
  if (!is.null(fea_tcell) && !is.null(fea_tcell$result) && nrow(fea_tcell$result) > 0) {
    write.csv(sanitize_gost_result(fea_tcell$result), file.path(results_dir, "fea_tcell.csv"), row.names = FALSE)
    cat("## T-Cells (Cluster 2) Enrichment\n\n", file = summary_file, append = TRUE)
    top_tcell <- head(fea_tcell$result %>% select(term_id, term_name, p_value) %>% arrange(p_value), 10)
    cat("Top 10 enriched Gene Ontology (Biological Process) terms:\n\n", file = summary_file, append = TRUE)
    cat(knitr::kable(top_tcell), sep = "\n", file = summary_file, append = TRUE)
    cat("\n\n", file = summary_file, append = TRUE)
  } else {
    cat("## T-Cells (Cluster 2) Enrichment\n\nNo enrichment results available (check log for details).\n\n", 
        file = summary_file, append = TRUE)
  }
  
  message("=== Assignment Plots and Tables Construction Completed ===")
}
