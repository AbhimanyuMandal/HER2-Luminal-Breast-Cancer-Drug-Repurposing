# ============================================================
# GSE176078 — Cluster Marker Validation
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

input_file <- "data/processed/GSE176078_RPCA_UMAP_clustered.rds"

output_dir <- "results/markers"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load object
# ------------------------------------------------------------

message("==============================================")
message("Loading clustered object...")
message("==============================================")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))
message("Clusters: ", length(unique(Idents(obj))))

# ------------------------------------------------------------
# 2. Use RNA assay
# ------------------------------------------------------------

DefaultAssay(obj) <- "RNA"

message("")
message("Default assay: ", DefaultAssay(obj))

# ------------------------------------------------------------
# 3. Join Seurat v5 RNA layers
# ------------------------------------------------------------

message("")
message("==============================================")
message("Joining RNA layers...")
message("==============================================")

obj <- JoinLayers(
  object = obj,
  assay = "RNA"
)

message("RNA layers joined successfully.")

# ------------------------------------------------------------
# 4. Cluster sizes
# ------------------------------------------------------------

message("")
message("==============================================")
message("Cluster sizes")
message("==============================================")

cluster_sizes <- as.data.frame(
  table(Idents(obj))
)

colnames(cluster_sizes) <- c(
  "cluster",
  "cells"
)

cluster_sizes <- cluster_sizes %>%
  arrange(desc(cells))

print(cluster_sizes)

write.csv(
  cluster_sizes,
  file.path(
    output_dir,
    "cluster_sizes.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Cluster × major cell type
# ------------------------------------------------------------

message("")
message("Generating cluster × major cell type table...")

if ("celltype_major" %in% colnames(obj@meta.data)) {

  cluster_celltype <- table(
    obj$seurat_clusters,
    obj$celltype_major
  )

  write.csv(
    as.data.frame(cluster_celltype),
    file.path(
      output_dir,
      "cluster_by_major_celltype.csv"
    ),
    row.names = FALSE
  )

} else {

  message("celltype_major not found in metadata.")

}

# ------------------------------------------------------------
# 6. Cluster × breast cancer subtype
# ------------------------------------------------------------

if ("subtype" %in% colnames(obj@meta.data)) {

  message("Generating cluster × subtype table...")

  cluster_subtype <- table(
    obj$seurat_clusters,
    obj$subtype
  )

  write.csv(
    as.data.frame(cluster_subtype),
    file.path(
      output_dir,
      "cluster_by_subtype.csv"
    ),
    row.names = FALSE
  )

} else {

  message("subtype not found in metadata.")

}

# ------------------------------------------------------------
# 7. Find cluster markers
# ------------------------------------------------------------

message("")
message("==============================================")
message("Finding cluster markers...")
message("==============================================")

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

# ------------------------------------------------------------
# 8. Check marker result
# ------------------------------------------------------------

if (nrow(markers) == 0) {

  stop(
    "No marker genes were identified. ",
    "Check RNA layers and normalization."
  )

}

message("")
message(
  "Marker genes identified: ",
  nrow(markers)
)

# ------------------------------------------------------------
# 9. Sort markers
# ------------------------------------------------------------

markers <- markers %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

write.csv(
  markers,
  file.path(
    output_dir,
    "all_cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Top 20 markers per cluster
# ------------------------------------------------------------

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top_markers,
  file.path(
    output_dir,
    "top20_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Print top markers
# ------------------------------------------------------------

message("")
message("==============================================")
message("TOP MARKERS")
message("==============================================")

for (cl in sort(unique(top_markers$cluster))) {

  message("")
  message("CLUSTER ", cl)

  print(
    top_markers %>%
      filter(cluster == cl) %>%
      select(
        gene,
        avg_log2FC,
        pct.1,
        pct.2,
        p_val_adj
      ) %>%
      head(10)
  )
}

# ------------------------------------------------------------
# 12. Save joined object
# ------------------------------------------------------------

joined_file <- "data/processed/GSE176078_RPCA_UMAP_clustered_joined.rds"

saveRDS(
  obj,
  joined_file,
  compress = FALSE
)

message("")
message("==============================================")
message("MARKER ANALYSIS COMPLETE")
message("==============================================")

message(
  "Markers: ",
  file.path(
    output_dir,
    "all_cluster_markers.csv"
  )
)

message(
  "Top markers: ",
  file.path(
    output_dir,
    "top20_markers_per_cluster.csv"
  )
)

message(
  "Joined object: ",
  joined_file
)