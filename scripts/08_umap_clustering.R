# ============================================================
# GSE176078 — UMAP and Clustering after Seurat RPCA
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

options(future.globals.maxSize = 6 * 1024^3)

input_file <- "data/processed/GSE176078_RPCA_integrated.rds"

output_file <- "data/processed/GSE176078_RPCA_UMAP_clustered.rds"

results_dir <- "results/umap"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

message("Loading RPCA-integrated object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))

# ------------------------------------------------------------
# Verify integrated reduction
# ------------------------------------------------------------

if (!"integrated.rpca" %in% Reductions(obj)) {
  stop("integrated.rpca reduction not found.")
}

# ------------------------------------------------------------
# Find neighbors
# ------------------------------------------------------------

message("")
message("Finding neighbors using integrated RPCA...")

obj <- FindNeighbors(
  obj,
  reduction = "integrated.rpca",
  dims = 1:30,
  k.param = 20,
  verbose = TRUE
)

# ------------------------------------------------------------
# Clustering
# ------------------------------------------------------------

message("")
message("Clustering cells...")

obj <- FindClusters(
  obj,
  resolution = 0.5,
  algorithm = 1,
  verbose = TRUE
)

message("")
message("Number of clusters: ",
        length(unique(Idents(obj))))

# ------------------------------------------------------------
# UMAP
# ------------------------------------------------------------

message("")
message("Running UMAP...")

obj <- RunUMAP(
  obj,
  reduction = "integrated.rpca",
  dims = 1:30,
  reduction.name = "umap.rpca",
  reduction.key = "UMAP_RPCA_",
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  verbose = TRUE
)

# ------------------------------------------------------------
# Plot: clusters
# ------------------------------------------------------------

message("Generating cluster UMAP...")

p_cluster <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("GSE176078 — Seurat RPCA Clusters")

ggsave(
  file.path(
    results_dir,
    "UMAP_RPCA_clusters.png"
  ),
  p_cluster,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Plot: major cell type
# ------------------------------------------------------------

message("Generating cell-type UMAP...")

p_celltype <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = "celltype_major",
  raster = TRUE
) +
  ggtitle("GSE176078 — Major Cell Types")

ggsave(
  file.path(
    results_dir,
    "UMAP_RPCA_celltypes.png"
  ),
  p_celltype,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Plot: patient
# ------------------------------------------------------------

message("Generating patient UMAP...")

p_patient <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = "patient_id",
  raster = TRUE
) +
  ggtitle("GSE176078 — Patient Mixing")

ggsave(
  file.path(
    results_dir,
    "UMAP_RPCA_patients.png"
  ),
  p_patient,
  width = 12,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# Plot: subtype
# ------------------------------------------------------------

message("Generating subtype UMAP...")

p_subtype <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = "subtype",
  raster = TRUE
) +
  ggtitle("GSE176078 — Breast Cancer Subtype")

ggsave(
  file.path(
    results_dir,
    "UMAP_RPCA_subtypes.png"
  ),
  p_subtype,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Cluster counts
# ------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(obj$seurat_clusters)
)

colnames(cluster_counts) <- c(
  "cluster",
  "cells"
)

write.csv(
  cluster_counts,
  file.path(
    results_dir,
    "cluster_cell_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Save object
# ------------------------------------------------------------

message("")
message("Saving clustered object...")

saveRDS(
  obj,
  output_file,
  compress = FALSE
)

message("")
message("==============================================")
message("UMAP + CLUSTERING COMPLETE")
message("==============================================")
message("Clusters: ", length(unique(obj$seurat_clusters)))
message("Saved: ", output_file)