# Breast Cancer Single-Cell RNA-seq (scRNA-seq) Pipeline

This dockerized project executes a modular single-cell RNA sequencing analysis pipeline for breast cancer tissue data. It starts from a pre-processed and normalized Seurat dataset, implements advanced downstream analysis steps (including clustering at multiple resolutions, marker detection, functional enrichment analysis, and cell cycle scoring), and generates publication-grade figures.

## Project Structure

```
breast-cancer-pipeline/
├── Dockerfile                  # Builds R environment with Seurat, gprofiler2, etc.
├── docker-compose.yml          # Container configuration and workspace mounting
├── main.R                      # Main execution script orchestrating the pipeline
├── README.md                   # This instruction and documentation file
├── data/                       # Directory containing the input datasets (ignored by git)
│   └── bc                      # Pre-processed breast cancer Seurat RDS object
├── src/                        # Modular analysis scripts
│   ├── data_loader.R           # Verifying data paths and reading RDS objects
│   ├── processing.R            # Clustering, marker detection, FEA, and cell cycle scoring
│   └── visualizer.R            # Exporting plots and FEA result tables
└── results/                    # Output directory for plots and final RDS object (mounted to host)
```

## Analytical Processes Performed

1. **Dimensionality Reduction**: Computes a 2D UMAP embedding based on the pre-computed PCA coordinates.
2. **Clustering**: Runs Seurat's graph-based clustering at three resolutions:
   - **0.3**: Coarser cell groupings.
   - **0.6**: Moderate clustering.
   - **1.2**: Finer, highly resolved subpopulation clusters.
3. **Marker Detection**: Performs differential expression (DE) analysis (`FindAllMarkers`) on the 0.6 resolution clusters to locate distinct gene markers for each group.
4. **Functional Enrichment Analysis (FEA)**: Performs enrichment checks via `gprofiler2` querying the `GO:BP` (Gene Ontology: Biological Process) database for:
   - **Cancer/Epithelial cells** (Cluster 12) markers.
   - **T-cells** (Cluster 2) markers.
5. **Cell Cycle Scoring**: Assigns S/G2M phase scores using Seurat's `CellCycleScoring` and identifies cells undergoing active proliferation.

## Generated Charts and Outputs (saved in `results/`)

All figures are automatically saved as high-resolution PNGs:

- **`umap_initial.png`**: Original UMAP projection of the cell communities.
- **`umap_clusters_res_0.6.png`**: UMAP showing the clustering at resolution 0.6.
- **`umap_clusters_res_0.3.png`** & **`umap_clusters_res_1.2.png`**: UMAPs showing alternative resolutions.
- **`umap_clusters_res_0.3_1.2.png`**: Multi-panel UMAP comparing resolution 0.3 and 1.2.
- **`umap_celltypes.png`**: UMAP projection labeled with annotated cell type names.
- **`featureplot_markers_MS4A1_KIT_CD14.png`**: Feature plot showing canonical markers (MS4A1 for B-cells, KIT for Mast cells, CD14 for Monocytes/Macrophages).
- **`vlnplot_markers_MS4A1_KIT_CD14.png`**: Violin plot of marker expressions grouped by annotated cell types.
- **`featureplot_top1_markers_per_cluster.png`**: Grid FeaturePlot showing the top marker gene for each cluster.
- **`heatmap_top6_markers.png`**: DoHeatmap showing the expression profiles of the top 6 markers for each cell type.
- **`vlnplot_cellcycle_scores.png`**: Violin plot showing S and G2M phase scores across all cell types.
- **`featureplot_cellcycle_scores.png`**: Comparative UMAP projection of cell cycle S.Score and G2M.Score.
- **`patchwork_q12.png`**: Combined diagnostic figure (CD3E feature plot, Cell Cycle Phase dim plot, and MKI67 violin plot).
- **`barchart_cellcycle_proportions.png`**: Dodge bar chart representing cell cycle phase percentages (G1, S, G2M) across cell types.
- **`fea_cancer_epithelial.csv`** & **`fea_tcell.csv`**: Raw enrichment table outputs.
- **`fea_summary.md`**: Markdown document summarizing the top 10 enriched terms for Cluster 12 and Cluster 2.
- **`processed_breast_cancer.rds`**: Serialized R object containing all added UMAP coordinates, clusters, and metadata.

---

## How to Get Up and Running

### Prerequisites
Make sure you have [Docker](https://www.docker.com/) and [Docker Compose](https://docs.docker.com/compose/) installed on your machine.

### Data Setup
Before executing the pipeline, you need to download and set up the pre-processed breast cancer dataset:
1. Download the data sample package from the [Breast Cancer scRNA-seq Dataset Link](https://1drv.ms/u/c/30fcedf48cd35cea/IQC2roPcIGrCQpMr7Z-wwFuxAZ0A59o0ExXrH0k_bfgs7pM?e=x1nJru).
2. Unzip the downloaded archive and place the dataset file (specifically `bc`) into the `data/` folder in the project root directory.

### Execution

To run the pipeline and generate the figures, execute the following commands in your terminal from the project root directory:

```bash
# Build the Docker image containing R, Seurat, and the dependencies
docker-compose build

# Start the container and run the pipeline
docker-compose up
```

Once the execution finishes, you will see a message:
`=== Pipeline Run Finished Successfully ===`

All generated figures, tables, and the processed Seurat object will be saved directly in the `results/` folder on your host machine.
