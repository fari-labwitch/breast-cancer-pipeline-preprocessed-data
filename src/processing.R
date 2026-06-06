# src/processing.R

library(Seurat)
library(dplyr)
library(gprofiler2)

#' Process Breast Cancer scRNA-seq Seurat Object
#'
#' @param se A Seurat object (already normalized and run with PCA)
#' @return A processed Seurat object containing UMAP, multiple clustering levels,
#'         differential markers, FEA results, and cell cycle scoring.
process_breast_cancer <- function(se) {
  message("=== Starting scRNA-seq Analysis Pipeline ===")
  
  # 1. Run UMAP embedding
  message("Running UMAP dimensionality reduction (dims 1:30)...")
  se <- RunUMAP(se, dims = 1:30)
  
  # 2. Run Clustering at multiple resolutions
  message("Finding neighbors in PCA space (dims 1:30)...")
  se <- FindNeighbors(se, dims = 1:30, reduction = "pca")
  
  message("Finding clusters at resolution 0.6...")
  se <- FindClusters(se, resolution = 0.6)
  
  message("Finding clusters at resolution 0.3...")
  se <- FindClusters(se, resolution = 0.3)
  
  message("Finding clusters at resolution 1.2...")
  se <- FindClusters(se, resolution = 1.2)
  
  # 3. Detect Marker Genes for Resolution 0.6
  message("Running marker detection for clusters at resolution 0.6...")
  se <- SetIdent(se, value = "SCT_snn_res.0.6")
  
  # Run FindAllMarkers to find cluster markers
  markers <- FindAllMarkers(
    se, 
    logfc.threshold = 0.5, 
    only.pos = TRUE, 
    max.cells.per.ident = 50, 
    verbose = FALSE
  )
  
  # Save markers inside Seurat object misc slot for visualization steps
  se@misc$markers <- markers
  
  # Print preview of top markers
  message("Top 2 markers per cluster by avg_log2FC:")
  top2_preview <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = 2, wt = avg_log2FC)
  print(head(top2_preview, 10))
  
  # 4. Functional Enrichment Analysis (FEA)
  # FEA for Cancer/Epithelial Cells (Cluster 12)
  message("Running FEA for Cancer/Epithelial cells (Cluster 12)...")
  cancer_markers <- markers %>%
    filter(cluster == "12", p_val_adj < 0.01, avg_log2FC > 1) %>%
    pull(gene)
  
  fea_cancer <- NULL
  if (length(cancer_markers) > 0) {
    message(sprintf("Found %d markers for Cluster 12. Querying g:Profiler...", length(cancer_markers)))
    fea_cancer <- tryCatch({
      gost(query = cancer_markers, 
           organism = "hsapiens", 
           significant = TRUE, 
           sources = "GO:BP")
    }, error = function(e) {
      warning("FEA query for Cluster 12 failed (check internet connection): ", e$message)
      NULL
    })
  } else {
    message("No markers for Cluster 12 met the criteria (p_val_adj < 0.01, avg_log2FC > 1).")
  }
  se@misc$fea_cancer <- fea_cancer
  
  # FEA for T-cells (Cluster 2)
  message("Running FEA for T-cells (Cluster 2)...")
  tcell_markers <- markers %>%
    filter(cluster == "2", avg_log2FC > 0.5) %>%
    pull(gene)
  
  fea_tcell <- NULL
  if (length(tcell_markers) > 0) {
    message(sprintf("Found %d markers for Cluster 2. Querying g:Profiler...", length(tcell_markers)))
    fea_tcell <- tryCatch({
      gost(query = tcell_markers, 
           organism = "hsapiens", 
           significant = TRUE, 
           sources = "GO:BP")
    }, error = function(e) {
      warning("FEA query for Cluster 2 failed (check internet connection): ", e$message)
      NULL
    })
  } else {
    message("No markers for Cluster 2 met the criteria (avg_log2FC > 0.5).")
  }
  se@misc$fea_tcell <- fea_tcell
  
  # 5. Cell Cycle Scoring
  message("Performing cell cycle scoring...")
  # cc.genes.updated.2019 is pre-loaded in Seurat
  se <- CellCycleScoring(
    se, 
    s.features = cc.genes.updated.2019$s.genes, 
    g2m.features = cc.genes.updated.2019$g2m.genes, 
    set.ident = TRUE
  )
  
  message("=== Analysis Pipeline Completed Successfully ===")
  return(se)
}
