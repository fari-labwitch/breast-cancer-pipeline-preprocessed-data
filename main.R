# main.R
source("src/data_loader.R")
source("src/processing.R")
source("src/visualizer.R")

message("=============================================")
message("=== Starting Breast Cancer scRNA Pipeline ===")
message("=============================================")

results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Run execution stages
raw_paths <- download_breast_cancer_data()
se_merged <- load_and_merge_samples(raw_paths)

se_processed <- process_breast_cancer(se_merged)

message("=== Commencing Figure Construction for Course Questions ===")
generate_assignment_plots(se_processed, results_dir = results_dir)

# Save intermediate Seurat object snapshot for later downstream parts (Part 8+)
message("Saving structured object snapshot...")
saveRDS(se_processed, file = file.path(results_dir, "processed_breast_cancer.rds"))

message("=== Pipeline Run Finished Successfully ===")