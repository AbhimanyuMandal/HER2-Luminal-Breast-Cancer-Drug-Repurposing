# ============================================================
# GSE113196 — NORMAL BREAST EPITHELIAL DATA
# Script 16: Normalization, HVG Selection & PCA
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

message("==============================================")
message("GSE113196 — NORMALIZATION + PCA")
message("==============================================")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/processed/GSE113196_normal_breast_QC_v1.rds"

output_dir <- "results/normal_breast"

output_file <-
  "data/processed/GSE113196_normal_breast_PCA_v1.rds"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load object
# ------------------------------------------------------------

message("")
message("Loading QC-filtered object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))

message("")
message("Individuals:")
print(table(obj$sample_id))

# ------------------------------------------------------------
# 3. Normalize RNA expression
# ------------------------------------------------------------

message("")
message("==============================================")
message("NORMALIZATION")
message("==============================================")

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

message("Normalization complete.")

# ------------------------------------------------------------
# 4. Identify highly variable genes
# ------------------------------------------------------------

message("")
message("==============================================")
message("HIGHLY VARIABLE GENES")
message("==============================================")

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

hvg <- VariableFeatures(obj)

message("Number of HVGs: ", length(hvg))

message("")
message("Top 30 HVGs:")

print(head(hvg, 30))

# Save HVGs

hvg_table <- data.frame(
  rank = seq_along(hvg),
  gene = hvg
)

write.csv(
  hvg_table,
  file.path(
    output_dir,
    "normal_breast_HVGs.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Scale data
# ------------------------------------------------------------

message("")
message("==============================================")
message("SCALING")
message("==============================================")

obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = TRUE
)

message("Scaling complete.")

# ------------------------------------------------------------
# 6. PCA
# ------------------------------------------------------------

message("")
message("==============================================")
message("PCA")
message("==============================================")

obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 50,
  verbose = TRUE
)

message("PCA complete.")

# ------------------------------------------------------------
# 7. PCA variance / elbow plot
# ------------------------------------------------------------

message("")
message("Generating ElbowPlot...")

p_elbow <- ElbowPlot(
  obj,
  ndims = 50
) +
  ggtitle(
    "GSE113196 Normal Breast — PCA Elbow Plot"
  )

ggsave(
  file.path(
    output_dir,
    "normal_breast_PCA_ElbowPlot.png"
  ),
  p_elbow,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 8. PCA loadings
# ------------------------------------------------------------

message("")
message("Top PCA-associated genes:")

for (pc in 1:10) {

  message("")
  message("PC ", pc)

loading_df <- as.data.frame(
  Loadings(obj)[, pc, drop = FALSE]
)

colnames(loading_df) <- "loading"

loading_df$gene <- rownames(loading_df)

loading_df <- loading_df %>%
  select(gene, loading) %>%
  mutate(
    absolute_loading = abs(loading)
  ) %>%
  arrange(
    desc(absolute_loading)
  )

print(
  head(
    loading_df,
    15
  )
)
}

# ------------------------------------------------------------
# 9. PCA dimensions summary
# ------------------------------------------------------------

pca_summary <- data.frame(
  PC = paste0("PC_", 1:50),
  stdev = Stdev(obj, reduction = "pca")[1:50]
)

pca_summary$variance_percent <-
  100 *
  pca_summary$stdev^2 /
  sum(pca_summary$stdev^2)

write.csv(
  pca_summary,
  file.path(
    output_dir,
    "normal_breast_PCA_variance.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Save object
# ------------------------------------------------------------

message("")
message("==============================================")
message("SAVING PCA OBJECT")
message("==============================================")

saveRDS(
  obj,
  output_file,
  compress = FALSE
)

message("")
message("Saved:")
message(output_file)

message("")
message("Additional outputs:")
message(
  file.path(
    output_dir,
    "normal_breast_HVGs.csv"
  )
)

message(
  file.path(
    output_dir,
    "normal_breast_PCA_ElbowPlot.png"
  )
)

message(
  file.path(
    output_dir,
    "normal_breast_PCA_variance.csv"
  )
)

message("")
message("==============================================")
message("SCRIPT 16 COMPLETE")
message("==============================================")