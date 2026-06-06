# src/data_loader.R

#' Download/Verify Breast Cancer scRNA-seq Data
#'
#' @return Character path to the dataset file
download_breast_cancer_data <- function() {
  dataset_path <- "data/bc"
  
  message("Verifying dataset path...")
  if (!file.exists(dataset_path)) {
    stop(paste("Dataset file not found at:", dataset_path, 
               "\nPlease ensure data/bc is placed in the project root directory."))
  }
  
  file_info <- file.info(dataset_path)
  size_mb <- round(file_info$size / (1024^2), 2)
  message(sprintf("Dataset verified: %s (%s MB)", dataset_path, size_mb))
  
  return(dataset_path)
}

#' Load the Seurat Object
#'
#' @param file_path Character path to the RDS file
#' @return A Seurat object
load_and_merge_samples <- function(file_path) {
  message(sprintf("Loading Seurat object from %s ...", file_path))
  bc <- readRDS(file_path)
  message("Seurat object loaded successfully.")
  print(bc)
  return(bc)
}
