library(Seurat)
library(dplyr)
library(ggplot2)

message("==============================================")
message("GSE113196 — INTEGRATED CLUSTERING")
message("==============================================")

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

input_file <- "data/processed/GSE113196_normal_breast_RPCA_v1.rds"

output_file <- "data/processed/GSE113196_normal_breast_RPCA_UMAP_clustered_v1.rds"

output_dir <- "results/normal_breast"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---------------------------------------------------------
# 1. Load integrated object
# ---------------------------------------------------------

message("")
message("Loading RPCA-integrated object...")

x <- readRDS(input_file)

message("Cells: ", ncol(x))
message("Genes: ", nrow(x))

message("")
message("Individuals:")
print(table(x$sample_id))

# ---------------------------------------------------------
# 2. Confirm PCA
# ---------------------------------------------------------

if (!"pca" %in% Reductions(x)) {
  stop("PCA reduction not found.")
}

message("")
message("Using PCA dimensions: 1-30")

# ---------------------------------------------------------
# 3. Find neighbors using integrated PCA
# ---------------------------------------------------------

message("")
message("Finding nearest neighbors...")

x <- FindNeighbors(
  x,
  reduction = "pca",
  dims = 1:30,
  graph.name = "integrated_snn",
  verbose = TRUE
)

message("Neighbor graph complete.")

# ---------------------------------------------------------
# 4. Clustering
# ---------------------------------------------------------

message("")
message("Finding clusters...")

x <- FindClusters(
  x,
  graph.name = "integrated_snn",
  resolution = 0.4,
  algorithm = 1,
  random.seed = 1234,
  verbose = TRUE
)

message("")
message(
  "Clusters identified: ",
  length(levels(Idents(x)))
)

# ---------------------------------------------------------
# 5. Cluster sizes
# ---------------------------------------------------------

cluster_sizes <- as.data.frame(
  table(Idents(x))
)

colnames(cluster_sizes) <- c(
  "cluster",
  "cells"
)

cluster_sizes <- cluster_sizes %>%
  mutate(
    percentage = cells / sum(cells) * 100
  )

message("")
message("Cluster sizes:")
print(cluster_sizes)

write.csv(
  cluster_sizes,
  file.path(
    output_dir,
    "integrated_cluster_sizes.csv"
  ),
  row.names = FALSE
)

# ---------------------------------------------------------
# 6. UMAP by cluster
# ---------------------------------------------------------

message("")
message("Generating cluster UMAP...")

p_cluster <- DimPlot(
  x,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "GSE113196 — Normal Breast Integrated Clusters"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "UMAP_integrated_clusters.png"
  ),
  p_cluster,
  width = 11,
  height = 8,
  dpi = 300
)

# ---------------------------------------------------------
# 7. UMAP by individual
# ---------------------------------------------------------

message("")
message("Generating donor UMAP...")

p_individual <- DimPlot(
  x,
  reduction = "umap",
  group.by = "sample_id"
) +
  ggtitle(
    "GSE113196 — Integrated UMAP by Individual"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "UMAP_integrated_by_individual.png"
  ),
  p_individual,
  width = 11,
  height = 8,
  dpi = 300
)

# ---------------------------------------------------------
# 8. Cluster × individual composition
# ---------------------------------------------------------

message("")
message("Calculating cluster × individual composition...")

cluster_individual_counts <- table(
  x$seurat_clusters,
  x$sample_id
)

cluster_individual_percent <- prop.table(
  cluster_individual_counts,
  margin = 1
) * 100

message("")
message("Cluster × individual percentage:")
print(round(cluster_individual_percent, 1))

write.csv(
  as.data.frame.matrix(cluster_individual_counts),
  file.path(
    output_dir,
    "integrated_cluster_by_individual_counts.csv"
  )
)

write.csv(
  as.data.frame.matrix(round(cluster_individual_percent, 2)),
  file.path(
    output_dir,
    "integrated_cluster_by_individual_percent.csv"
  )
)

# ---------------------------------------------------------
# 9. Marker preparation
# ---------------------------------------------------------

message("")
message("Setting RNA assay for biological marker analysis...")

if ("RNA" %in% Assays(x)) {
  DefaultAssay(x) <- "RNA"
}

# ---------------------------------------------------------
# 10. Save
# ---------------------------------------------------------

message("")
message("Saving clustered integrated object...")

saveRDS(
  x,
  output_file
)

# ---------------------------------------------------------
# Final summary
# ---------------------------------------------------------

message("")
message("==============================================")
message("INTEGRATED CLUSTERING COMPLETE")
message("==============================================")

message("")
message("Cells: ", ncol(x))
message("Genes: ", nrow(x))
message(
  "Clusters: ",
  length(levels(Idents(x)))
)

message("")
message("Files saved:")
message(output_file)
message(
  file.path(
    output_dir,
    "integrated_cluster_sizes.csv"
  )
)
message(
  file.path(
    output_dir,
    "UMAP_integrated_clusters.png"
  )
)
message(
  file.path(
    output_dir,
    "UMAP_integrated_by_individual.png"
  )
)
message(
  file.path(
    output_dir,
    "integrated_cluster_by_individual_counts.csv"
  )
)
message(
  file.path(
    output_dir,
    "integrated_cluster_by_individual_percent.csv"
  )
)

message("")
message("==============================================")