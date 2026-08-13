# ============================================================
# GSE113196 — NORMAL BREAST MARKER IDENTIFICATION
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

message("==============================================")
message("GSE113196 — NORMAL BREAST MARKER ANALYSIS")
message("==============================================")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

input_file <-
  "data/processed/GSE113196_normal_breast_RPCA_UMAP_clustered_v1.rds"

output_dir <- "results/normal_breast/markers"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_object <-
  "data/processed/GSE113196_normal_breast_RPCA_UMAP_clustered_v1_markers.rds"

# ------------------------------------------------------------
# 1. Load object
# ------------------------------------------------------------

message("")
message("Loading clustered normal breast object...")

x <- readRDS(input_file)

message("Cells: ", ncol(x))
message("Genes: ", nrow(x))
message("Clusters: ", length(levels(Idents(x))))

# ------------------------------------------------------------
# 2. Set RNA assay
# ------------------------------------------------------------

if (!"RNA" %in% Assays(x)) {
  stop("RNA assay not found.")
}

DefaultAssay(x) <- "RNA"

message("")
message("Default assay: ", DefaultAssay(x))

# ------------------------------------------------------------
# 3. Join RNA layers
# ------------------------------------------------------------

message("")
message("Checking RNA assay layers...")

print(Layers(x[["RNA"]]))

message("")
message("Joining RNA layers for differential expression...")

x[["RNA"]] <- JoinLayers(x[["RNA"]])

message("RNA layers after JoinLayers:")

print(Layers(x[["RNA"]]))

# ------------------------------------------------------------
# 4. Confirm clusters
# ------------------------------------------------------------

message("")
message("Cluster sizes:")

cluster_sizes <- as.data.frame(
  table(Idents(x))
)

colnames(cluster_sizes) <- c(
  "cluster",
  "cells"
)

cluster_sizes$percentage <-
  cluster_sizes$cells /
  sum(cluster_sizes$cells) * 100

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
# 5. Find markers
# ------------------------------------------------------------

message("")
message("Finding cluster markers...")

markers <- FindAllMarkers(
  x,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

message("")
message("Markers identified: ", nrow(markers))

if (nrow(markers) == 0) {
  stop(
    paste(
      "No markers were identified.",
      "Check RNA layers and expression data."
    )
  )
}

# ------------------------------------------------------------
# 6. Save complete marker table
# ------------------------------------------------------------

write.csv(
  markers,
  file.path(
    output_dir,
    "all_cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Top 20 markers
# ------------------------------------------------------------

top20 <- markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top20,
  file.path(
    output_dir,
    "top20_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Top 10 markers
# ------------------------------------------------------------

top10 <- markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top10,
  file.path(
    output_dir,
    "top10_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Print top markers
# ------------------------------------------------------------

message("")
message("==============================================")
message("TOP MARKERS PER CLUSTER")
message("==============================================")

for (cl in levels(Idents(x))) {

  message("")
  message("Cluster ", cl)

  tmp <- markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 15)

  if (nrow(tmp) > 0) {

    print(
      tmp %>%
        select(
          gene,
          avg_log2FC,
          pct.1,
          pct.2,
          p_val_adj
        )
    )

  } else {

    message("No significant positive markers found.")

  }
}

# ------------------------------------------------------------
# 10. Canonical marker panel
# ------------------------------------------------------------

canonical_markers <- c(

  # General epithelial
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",

  # Luminal
  "KRT7",
  "MUC1",
  "KRT23",
  "KRT24",

  # Hormone receptor-associated
  "ESR1",
  "PGR",
  "AGR2",
  "GATA3",
  "FOXA1",

  # Luminal progenitor / secretory
  "KIT",
  "ALDH1A3",
  "KRT15",
  "KRT16",
  "KRT6A",
  "KRT6B",
  "S100A7",
  "SLPI",

  # Basal / myoepithelial
  "KRT5",
  "KRT14",
  "KRT17",
  "TP63",
  "ACTA2",
  "TAGLN",
  "MYL9",
  "CNN1",

  # Secretory
  "LALBA",
  "CSN1S1",
  "CSN2",
  "CSN3",
  "PRLR",

  # Proliferation
  "MKI67",
  "TOP2A",
  "UBE2C",
  "BIRC5",
  "STMN1",

  # Stromal contamination
  "COL1A1",
  "COL1A2",
  "DCN",
  "LUM",

  # Endothelial contamination
  "PECAM1",
  "VWF",
  "EMCN",
  "KDR",

  # Immune contamination
  "PTPRC",
  "CD3D",
  "CD3E",
  "MS4A1",
  "CD79A",
  "LYZ"
)

canonical_markers <- intersect(
  canonical_markers,
  rownames(x)
)

message("")
message(
  "Canonical markers available: ",
  length(canonical_markers)
)

# ------------------------------------------------------------
# 11. Canonical marker DotPlot
# ------------------------------------------------------------

message("")
message("Generating canonical marker DotPlot...")

p_dot <- DotPlot(
  x,
  features = canonical_markers,
  assay = "RNA",
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle(
    "GSE113196 — Canonical Marker Expression"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "canonical_marker_DotPlot.png"
  ),
  p_dot,
  width = 16,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 12. Key FeaturePlots
# ------------------------------------------------------------

key_markers <- intersect(
  c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT5",
    "KRT14",
    "KRT17",
    "KRT6A",
    "KRT6B",
    "MUC1",
    "KRT7",
    "ESR1",
    "PGR",
    "GATA3",
    "FOXA1",
    "KIT",
    "ALDH1A3",
    "MKI67",
    "TOP2A",
    "COL1A1",
    "DCN",
    "PECAM1",
    "PTPRC"
  ),
  rownames(x)
)

message("")
message("Generating key FeaturePlots...")

p_feature <- FeaturePlot(
  x,
  features = key_markers,
  reduction = "umap",
  ncol = 4
)

ggsave(
  file.path(
    output_dir,
    "key_epithelial_FeaturePlots.png"
  ),
  p_feature,
  width = 16,
  height = 14,
  dpi = 300
)

# ------------------------------------------------------------
# 13. Store marker results
# ------------------------------------------------------------

x@misc$normal_breast_marker_analysis <- list(
  markers = markers,
  top20 = top20,
  top10 = top10
)

# ------------------------------------------------------------
# 14. Save object
# ------------------------------------------------------------

message("")
message("Saving marker-characterized object...")

saveRDS(
  x,
  output_object
)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

message("")
message("==============================================")
message("MARKER ANALYSIS COMPLETE")
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
message(output_object)
message(
  file.path(
    output_dir,
    "all_cluster_markers.csv"
  )
)
message(
  file.path(
    output_dir,
    "top20_markers_per_cluster.csv"
  )
)
message(
  file.path(
    output_dir,
    "top10_markers_per_cluster.csv"
  )
)
message(
  file.path(
    output_dir,
    "canonical_marker_DotPlot.png"
  )
)
message(
  file.path(
    output_dir,
    "key_epithelial_FeaturePlots.png"
  )
)

message("")
message("==============================================")