# ============================================================
# GSE176078 — Normalization, HVG Selection & PCA
# Seurat v5
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(12345)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/processed/GSE176078_QC_filtered.rds"

output_dir <- "data/processed"
results_dir <- "results/pca"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load QC-filtered object
# ------------------------------------------------------------

message("Loading QC-filtered object...")

obj <- readRDS(input_file)

message(
  "Loaded: ",
  ncol(obj),
  " cells × ",
  nrow(obj),
  " genes"
)

# ------------------------------------------------------------
# 3. Normalize RNA expression
# ------------------------------------------------------------

message("")
message("Running LogNormalize...")

obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

# ------------------------------------------------------------
# 4. Identify highly variable genes
# ------------------------------------------------------------

message("")
message("Finding highly variable genes...")

obj <- FindVariableFeatures(
  obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

message(
  "Number of variable features: ",
  length(VariableFeatures(obj))
)

message("")
message("Top 20 variable genes:")

print(
  head(
    VariableFeatures(obj),
    20
  )
)

# ------------------------------------------------------------
# 5. Scale variable genes
# ------------------------------------------------------------

message("")
message("Scaling variable genes...")

obj <- ScaleData(
  obj,
  assay = "RNA",
  features = VariableFeatures(obj),
  verbose = TRUE
)

# ------------------------------------------------------------
# 6. Principal Component Analysis
# ------------------------------------------------------------

message("")
message("Running PCA...")

obj <- RunPCA(
  obj,
  assay = "RNA",
  features = VariableFeatures(obj),
  npcs = 50,
  verbose = TRUE
)

# ------------------------------------------------------------
# 7. Save PCA object
# ------------------------------------------------------------

message("")
message("Saving object...")

saveRDS(
  obj,
  file.path(
    output_dir,
    "GSE176078_QC_normalized_PCA.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# 8. Elbow plot
# ------------------------------------------------------------

message("Generating elbow plot...")

p_elbow <- ElbowPlot(
  obj,
  ndims = 50
) +
  ggtitle(
    "GSE176078 — PCA Elbow Plot"
  )

ggsave(
  file.path(
    results_dir,
    "PCA_elbow_plot.png"
  ),
  p_elbow,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 9. PCA variance summary
# ------------------------------------------------------------

pca_stdev <- Stdev(
  obj[["pca"]]
)

pca_variance <- pca_stdev^2

pca_variance_percent <- (
  pca_variance /
    sum(pca_variance)
) * 100

pca_summary <- data.frame(
  PC = seq_along(pca_stdev),
  StandardDeviation = pca_stdev,
  Variance = pca_variance,
  VariancePercent = pca_variance_percent
)

write.csv(
  pca_summary,
  file.path(
    results_dir,
    "PCA_variance_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Final summary
# ------------------------------------------------------------

message("")
message("==========================================")
message("NORMALIZATION + PCA COMPLETE")
message("==========================================")

message(
  "Cells: ",
  ncol(obj)
)

message(
  "Genes: ",
  nrow(obj)
)

message(
  "HVGs: ",
  length(VariableFeatures(obj))
)

message(
  "PCA dimensions: ",
  ncol(
    Embeddings(obj, "pca")
  )
)

message("")
message(
  "Saved: ",
  file.path(
    output_dir,
    "GSE176078_QC_normalized_PCA.rds"
  )
)

message(
  "Elbow plot: ",
  file.path(
    results_dir,
    "PCA_elbow_plot.png"
  )
)

message("==========================================")