library(Seurat)
library(dplyr)

message("========================================")
message("ANNOTATION VALIDATION")
message("========================================")

input_file <- "data/processed/GSE176078_RPCA_UMAP_clustered_joined.rds"

message("Loading clustered object...")
obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Clusters: ", length(unique(Idents(obj))))

# ---------------------------------------------------------
# 1. Inspect metadata
# ---------------------------------------------------------

message("")
message("Metadata columns:")
print(colnames(obj@meta.data))

# ---------------------------------------------------------
# 2. Check original major cell-type annotations
# ---------------------------------------------------------

if ("celltype_major" %in% colnames(obj@meta.data)) {

  message("")
  message("========================================")
  message("ORIGINAL CELLTYPE MAJOR")
  message("========================================")

  print(table(obj$celltype_major))

} else {

  message("")
  message("WARNING: celltype_major not found.")
}

# ---------------------------------------------------------
# 3. Cluster x original cell type
# ---------------------------------------------------------

if ("celltype_major" %in% colnames(obj@meta.data)) {

  message("")
  message("========================================")
  message("CLUSTER × ORIGINAL CELL TYPE")
  message("========================================")

  cluster_celltype <- table(
    Cluster = Idents(obj),
    CellType = obj$celltype_major
  )

  print(cluster_celltype)

  # Save counts
  write.csv(
    as.data.frame(cluster_celltype),
    "results/markers/cluster_vs_original_celltype_counts.csv",
    row.names = FALSE
  )

  # Row percentages
  cluster_pct <- prop.table(cluster_celltype, margin = 1) * 100

  message("")
  message("========================================")
  message("CLUSTER × ORIGINAL CELL TYPE (%)")
  message("========================================")

  print(round(cluster_pct, 1))

  write.csv(
    as.data.frame(cluster_pct),
    "results/markers/cluster_vs_original_celltype_percent.csv",
    row.names = FALSE
  )
}

# ---------------------------------------------------------
# 4. Cluster sizes
# ---------------------------------------------------------

message("")
message("========================================")
message("CLUSTER SIZES")
message("========================================")

cluster_sizes <- as.data.frame(table(Idents(obj)))

colnames(cluster_sizes) <- c("cluster", "cells")

print(cluster_sizes)

write.csv(
  cluster_sizes,
  "results/markers/cluster_sizes.csv",
  row.names = FALSE
)

# ---------------------------------------------------------
# 5. Print dominant original annotation per cluster
# ---------------------------------------------------------

if ("celltype_major" %in% colnames(obj@meta.data)) {

  message("")
  message("========================================")
  message("DOMINANT ORIGINAL ANNOTATION PER CLUSTER")
  message("========================================")

  dominant_annotation <- obj@meta.data %>%
    mutate(cluster = as.character(Idents(obj))) %>%
    count(cluster, celltype_major, name = "n") %>%
    group_by(cluster) %>%
    mutate(
      total = sum(n),
      percent = 100 * n / total
    ) %>%
    arrange(cluster, desc(n)) %>%
    slice_head(n = 3) %>%
    ungroup()

  print(dominant_annotation, n = Inf)

  write.csv(
    dominant_annotation,
    "results/markers/dominant_original_annotations.csv",
    row.names = FALSE
  )
}

# ---------------------------------------------------------
# 6. Final message
# ---------------------------------------------------------

message("")
message("========================================")
message("VALIDATION COMPLETE")
message("========================================")

message("Saved:")
message("results/markers/cluster_vs_original_celltype_counts.csv")
message("results/markers/cluster_vs_original_celltype_percent.csv")
message("results/markers/cluster_sizes.csv")
message("results/markers/dominant_original_annotations.csv")