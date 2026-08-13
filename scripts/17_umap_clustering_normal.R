library(Seurat)
library(dplyr)
library(ggplot2)

message("==============================================")
message("GSE113196 — NORMAL BREAST UMAP + CLUSTERING")
message("==============================================")

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

input_file <- "data/processed/GSE113196_normal_breast_PCA_v1.rds"

output_file <- "data/processed/GSE113196_normal_breast_UMAP_clustered_v1.rds"

output_dir <- "results/normal_breast"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 1. Load PCA object
# ---------------------------------------------------------

message("")
message("Loading PCA object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))

if ("sample_id" %in% colnames(obj@meta.data)) {
  
  message("")
  message("Individuals:")
  print(table(obj$sample_id))
  
}

# ---------------------------------------------------------
# 2. PCA dimensions
# ---------------------------------------------------------

dims_use <- 1:30

message("")
message("Using PCA dimensions: 1–30")
message("Number of PCs: ", length(dims_use))

# ---------------------------------------------------------
# 3. Find neighbors
# ---------------------------------------------------------

message("")
message("Finding nearest neighbors...")

obj <- FindNeighbors(
  obj,
  dims = dims_use,
  verbose = TRUE
)

message("Neighbor graph complete.")

# ---------------------------------------------------------
# 4. Clustering
# ---------------------------------------------------------

message("")
message("Finding clusters...")

obj <- FindClusters(
  obj,
  resolution = 0.5,
  verbose = TRUE
)

message("")
message("Clusters identified: ",
        length(unique(Idents(obj))))

message("")
message("Cluster sizes:")

cluster_sizes <- as.data.frame(table(Idents(obj)))

colnames(cluster_sizes) <- c(
  "cluster",
  "cells"
)

cluster_sizes <- cluster_sizes %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  )

print(cluster_sizes)

write.csv(
  cluster_sizes,
  file.path(
    output_dir,
    "cluster_sizes.csv"
  ),
  row.names = FALSE
)

# ---------------------------------------------------------
# 5. Generate UMAP
# ---------------------------------------------------------

message("")
message("Generating UMAP...")

set.seed(1234)

obj <- RunUMAP(
  obj,
  dims = dims_use,
  reduction = "pca",
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE
)

message("UMAP complete.")

# ---------------------------------------------------------
# 6. UMAP — clusters
# ---------------------------------------------------------

message("")
message("Generating cluster UMAP...")

p_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "GSE113196 — Normal Breast Epithelial Cells"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_clusters.png"
  ),
  plot = p_cluster,
  width = 10,
  height = 8,
  dpi = 300
)

# ---------------------------------------------------------
# 7. UMAP — individuals
# ---------------------------------------------------------

if ("sample_id" %in% colnames(obj@meta.data)) {
  
  message("")
  message("Generating individual/sample UMAP...")
  
  p_sample <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "sample_id"
  ) +
    ggtitle(
      "GSE113196 — UMAP by Individual"
    ) +
    theme_classic()
  
  ggsave(
    filename = file.path(
      output_dir,
      "UMAP_by_individual.png"
    ),
    plot = p_sample,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ---------------------------------------------------------
# 8. Cluster × individual composition
# ---------------------------------------------------------

if ("sample_id" %in% colnames(obj@meta.data)) {
  
  message("")
  message("Calculating cluster × individual composition...")
  
  cluster_sample <- table(
    Cluster = Idents(obj),
    Individual = obj$sample_id
  )
  
  write.csv(
    as.data.frame(cluster_sample),
    file.path(
      output_dir,
      "cluster_by_individual_counts.csv"
    ),
    row.names = FALSE
  )
  
  cluster_sample_pct <- prop.table(
    cluster_sample,
    margin = 1
  ) * 100
  
  write.csv(
    as.data.frame(cluster_sample_pct),
    file.path(
      output_dir,
      "cluster_by_individual_percent.csv"
    ),
    row.names = FALSE
  )
  
  message("")
  message("Cluster × individual composition:")
  print(round(cluster_sample_pct, 1))
}

# ---------------------------------------------------------
# 9. Save clustered object
# ---------------------------------------------------------

message("")
message("Saving clustered object...")

saveRDS(
  obj,
  output_file
)

# ---------------------------------------------------------
# 10. Final summary
# ---------------------------------------------------------

message("")
message("==============================================")
message("UMAP + CLUSTERING COMPLETE")
message("==============================================")

message("")
message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))
message(
  "Clusters: ",
  length(unique(Idents(obj)))
)

message("")
message("Files saved:")

message(output_file)

message(
  file.path(
    output_dir,
    "cluster_sizes.csv"
  )
)

message(
  file.path(
    output_dir,
    "UMAP_clusters.png"
  )
)

if ("sample_id" %in% colnames(obj@meta.data)) {
  
  message(
    file.path(
      output_dir,
      "UMAP_by_individual.png"
    )
  )
  
  message(
    file.path(
      output_dir,
      "cluster_by_individual_counts.csv"
    )
  )
  
  message(
    file.path(
      output_dir,
      "cluster_by_individual_percent.csv"
    )
  )
}

message("")
message("==============================================")