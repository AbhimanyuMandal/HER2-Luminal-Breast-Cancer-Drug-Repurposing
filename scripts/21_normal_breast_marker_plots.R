# ============================================================
# GSE113196 — NORMAL BREAST MARKER ANALYSIS + VISUALIZATION
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

input_file <-
  "data/processed/GSE113196_normal_breast_RPCA_UMAP_clustered_v1.rds"

output_dir <- "results/normal_breast/markers"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("==============================================")
message("GSE113196 — NORMAL BREAST MARKER ANALYSIS")
message("==============================================")

# ------------------------------------------------------------
# LOAD OBJECT
# ------------------------------------------------------------

x <- readRDS(input_file)

DefaultAssay(x) <- "RNA"

message("Cells: ", ncol(x))
message("Genes: ", nrow(x))
message("Clusters: ", length(levels(Idents(x))))

# ------------------------------------------------------------
# JOIN RNA LAYERS
# ------------------------------------------------------------

message("")
message("Joining RNA layers...")

x[["RNA"]] <- JoinLayers(x[["RNA"]])

message("RNA layers joined.")

# ------------------------------------------------------------
# FIND MARKERS
# ------------------------------------------------------------

message("")
message("Finding cluster markers...")

markers <- FindAllMarkers(
  x,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
)

message("")
message("Markers identified: ", nrow(markers))

# ------------------------------------------------------------
# SAVE MARKERS IMMEDIATELY
# ------------------------------------------------------------

write.csv(
  markers,
  file.path(
    output_dir,
    "all_cluster_markers.csv"
  ),
  row.names = FALSE
)

saveRDS(
  markers,
  file.path(
    output_dir,
    "all_cluster_markers.rds"
  )
)

message("Marker table saved.")

# ------------------------------------------------------------
# TOP MARKERS
# ------------------------------------------------------------

top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top_markers,
  file.path(
    output_dir,
    "top20_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

message("Top marker table saved.")

# ------------------------------------------------------------
# CANONICAL MARKER PANEL
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
  "ESR1",
  "PGR",
  "GATA3",
  "FOXA1",
  "AGR2",

  # Basal
  "KRT5",
  "KRT14",
  "KRT17",
  "KRT6A",
  "KRT6B",
  "TP63",

  # Secretory / epithelial states
  "ALDH1A3",
  "SCGB2A2",
  "LTF",
  "AQP5",
  "PIP",
  "BPIFB1",
  "SPDEF",
  "TFF3",

  # Cycling
  "MKI67",
  "TOP2A",
  "UBE2C",
  "BIRC5",
  "AURKB",

  # Stromal
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "DCN",
  "LUM",
  "SFRP4",
  "MFAP4",

  # Endothelial
  "PECAM1",
  "VWF",
  "EMCN",
  "KDR",

  # Immune
  "PTPRC",
  "CD3D",
  "CD3E",
  "MS4A1",
  "CD79A",
  "LYZ"
)

rna_genes <- rownames(x[["RNA"]])

canonical_available <- intersect(
  canonical_markers,
  rna_genes
)

message("")
message(
  "Canonical markers available: ",
  length(canonical_available)
)

print(canonical_available)

# ------------------------------------------------------------
# CANONICAL DOTPLOT
# ------------------------------------------------------------

if (length(canonical_available) > 0) {

  message("")
  message("Generating canonical marker DotPlot...")

  p_dot <- DotPlot(
    x,
    assay = "RNA",
    features = canonical_available,
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
    height = 10,
    dpi = 300
  )

  message("DotPlot saved.")

} else {

  message("WARNING: No canonical markers found.")

}

# ------------------------------------------------------------
# TOP MARKER HEATMAP
# ------------------------------------------------------------

heatmap_genes <- markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()

heatmap_genes <- intersect(
  heatmap_genes,
  rna_genes
)

message("")
message(
  "Generating heatmap with ",
  length(heatmap_genes),
  " genes..."
)

p_heatmap <- DoHeatmap(
  x,
  assay = "RNA",
  features = heatmap_genes,
  group.by = "seurat_clusters"
) +
  ggtitle(
    "GSE113196 — Top Cluster Markers"
  )

ggsave(
  file.path(
    output_dir,
    "top_cluster_markers_heatmap.png"
  ),
  p_heatmap,
  width = 16,
  height = 12,
  dpi = 300
)

message("Heatmap saved.")

# ------------------------------------------------------------
# SAVE UPDATED OBJECT
# ------------------------------------------------------------

saveRDS(
  x,
  file.path(
    output_dir,
    "GSE113196_normal_breast_marker_analysis_v1.rds"
  )
)

message("")
message("==============================================")
message("MARKER ANALYSIS COMPLETE")
message("==============================================")

message("")
message("Outputs:")
message("all_cluster_markers.csv")
message("all_cluster_markers.rds")
message("top20_markers_per_cluster.csv")
message("canonical_marker_DotPlot.png")
message("top_cluster_markers_heatmap.png")
message("GSE113196_normal_breast_marker_analysis_v1.rds")