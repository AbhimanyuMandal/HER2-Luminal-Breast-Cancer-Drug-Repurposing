# ============================================================
# GSE176078 — PCA Patient/Subtype Diagnostics
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

obj <- readRDS(
  "data/processed/GSE176078_QC_normalized_PCA.rds"
)

results_dir <- "results/pca"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# PCA colored by cancer subtype
# ------------------------------------------------------------

p_subtype <- DimPlot(
  obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "subtype",
  raster = TRUE
) +
  ggtitle("PCA — Cancer Subtype")

ggsave(
  file.path(
    results_dir,
    "PCA_PC1_PC2_by_subtype.png"
  ),
  p_subtype,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# PCA colored by major cell type
# ------------------------------------------------------------

p_celltype <- DimPlot(
  obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "celltype_major",
  raster = TRUE
) +
  ggtitle("PCA — Major Cell Type")

ggsave(
  file.path(
    results_dir,
    "PCA_PC1_PC2_by_celltype.png"
  ),
  p_celltype,
  width = 9,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# PCA colored by patient
# ------------------------------------------------------------

p_patient <- DimPlot(
  obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "patient_id",
  raster = TRUE
) +
  ggtitle("PCA — Patient")

ggsave(
  file.path(
    results_dir,
    "PCA_PC1_PC2_by_patient.png"
  ),
  p_patient,
  width = 10,
  height = 8,
  dpi = 300
)

message("")
message("PCA diagnostic plots saved.")