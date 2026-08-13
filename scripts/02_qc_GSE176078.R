# ============================================================
# GSE176078 — Quality Control Assessment
# 100,064 cells | 26 breast cancer patients
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(data.table)
})

set.seed(12345)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/processed/GSE176078_raw_seurat.rds"

qc_dir <- "results/qc"

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load Seurat object
# ------------------------------------------------------------

message("Loading Seurat object...")

obj <- readRDS(input_file)

message(
  "Loaded: ",
  ncol(obj),
  " cells × ",
  nrow(obj),
  " genes"
)

# ------------------------------------------------------------
# 3. QC metrics
# ------------------------------------------------------------

# Recalculate mitochondrial percentage independently
obj[["percent.mt.calc"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# Extract metadata only for QC plotting
qc <- obj@meta.data

qc$cell_id <- rownames(qc)

# ------------------------------------------------------------
# 4. Basic QC summary
# ------------------------------------------------------------

message("")
message("==========================================")
message("QC SUMMARY")
message("==========================================")

message("Cells: ", nrow(qc))

message("")
message("nFeature_RNA:")
print(summary(qc$nFeature_RNA))

message("")
message("nCount_RNA:")
print(summary(qc$nCount_RNA))

message("")
message("percent.mito:")
print(summary(qc$percent.mito))

message("")
message("percent.mt.calc:")
print(summary(qc$percent.mt.calc))

# ------------------------------------------------------------
# 5. Quantiles
# ------------------------------------------------------------

qc_quantiles <- data.frame(
  metric = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mito",
    "percent.mt.calc"
  ),
  Q01 = c(
    quantile(qc$nFeature_RNA, 0.01, na.rm = TRUE),
    quantile(qc$nCount_RNA, 0.01, na.rm = TRUE),
    quantile(qc$percent.mito, 0.01, na.rm = TRUE),
    quantile(qc$percent.mt.calc, 0.01, na.rm = TRUE)
  ),
  Q05 = c(
    quantile(qc$nFeature_RNA, 0.05, na.rm = TRUE),
    quantile(qc$nCount_RNA, 0.05, na.rm = TRUE),
    quantile(qc$percent.mito, 0.05, na.rm = TRUE),
    quantile(qc$percent.mt.calc, 0.05, na.rm = TRUE)
  ),
  Q25 = c(
    quantile(qc$nFeature_RNA, 0.25, na.rm = TRUE),
    quantile(qc$nCount_RNA, 0.25, na.rm = TRUE),
    quantile(qc$percent.mito, 0.25, na.rm = TRUE),
    quantile(qc$percent.mt.calc, 0.25, na.rm = TRUE)
  ),
  Median = c(
    median(qc$nFeature_RNA, na.rm = TRUE),
    median(qc$nCount_RNA, na.rm = TRUE),
    median(qc$percent.mito, na.rm = TRUE),
    median(qc$percent.mt.calc, na.rm = TRUE)
  ),
  Q75 = c(
    quantile(qc$nFeature_RNA, 0.75, na.rm = TRUE),
    quantile(qc$nCount_RNA, 0.75, na.rm = TRUE),
    quantile(qc$percent.mito, 0.75, na.rm = TRUE),
    quantile(qc$percent.mt.calc, 0.75, na.rm = TRUE)
  ),
  Q95 = c(
    quantile(qc$nFeature_RNA, 0.95, na.rm = TRUE),
    quantile(qc$nCount_RNA, 0.95, na.rm = TRUE),
    quantile(qc$percent.mito, 0.95, na.rm = TRUE),
    quantile(qc$percent.mt.calc, 0.95, na.rm = TRUE)
  ),
  Q99 = c(
    quantile(qc$nFeature_RNA, 0.99, na.rm = TRUE),
    quantile(qc$nCount_RNA, 0.99, na.rm = TRUE),
    quantile(qc$percent.mito, 0.99, na.rm = TRUE),
    quantile(qc$percent.mt.calc, 0.99, na.rm = TRUE)
  )
)

print(qc_quantiles)

fwrite(
  qc_quantiles,
  file.path(qc_dir, "qc_quantiles.csv")
)

# ------------------------------------------------------------
# 6. Patient-level QC
# ------------------------------------------------------------

patient_qc <- as.data.frame(
  aggregate(
    cbind(
      nFeature_RNA,
      nCount_RNA,
      percent.mito
    ) ~ patient_id,
    data = qc,
    FUN = median
  )
)

names(patient_qc) <- c(
  "patient_id",
  "median_nFeature_RNA",
  "median_nCount_RNA",
  "median_percent_mito"
)

patient_qc$cell_count <- as.integer(
  table(qc$patient_id)[patient_qc$patient_id]
)

patient_qc <- patient_qc[
  order(-patient_qc$cell_count),
]

fwrite(
  patient_qc,
  file.path(qc_dir, "patient_qc_summary.csv")
)

# ------------------------------------------------------------
# 7. Subtype-level QC
# ------------------------------------------------------------

subtype_qc <- as.data.frame(
  aggregate(
    cbind(
      nFeature_RNA,
      nCount_RNA,
      percent.mito
    ) ~ subtype,
    data = qc,
    FUN = median
  )
)

names(subtype_qc) <- c(
  "subtype",
  "median_nFeature_RNA",
  "median_nCount_RNA",
  "median_percent_mito"
)

subtype_qc$cell_count <- as.integer(
  table(qc$subtype)[subtype_qc$subtype]
)

fwrite(
  subtype_qc,
  file.path(qc_dir, "subtype_qc_summary.csv")
)

# ------------------------------------------------------------
# 8. QC plot — nFeature_RNA
# ------------------------------------------------------------

p1 <- ggplot(
  qc,
  aes(x = nFeature_RNA)
) +
  geom_histogram(
    bins = 100
  ) +
  scale_x_log10() +
  labs(
    title = "GSE176078 — Detected Genes per Cell",
    x = "nFeature_RNA",
    y = "Number of cells"
  ) +
  theme_classic()

ggsave(
  file.path(qc_dir, "qc_nFeature_distribution.png"),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 9. QC plot — nCount_RNA
# ------------------------------------------------------------

p2 <- ggplot(
  qc,
  aes(x = nCount_RNA)
) +
  geom_histogram(
    bins = 100
  ) +
  scale_x_log10() +
  labs(
    title = "GSE176078 — RNA Counts per Cell",
    x = "nCount_RNA",
    y = "Number of cells"
  ) +
  theme_classic()

ggsave(
  file.path(qc_dir, "qc_nCount_distribution.png"),
  p2,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 10. QC plot — mitochondrial percentage
# ------------------------------------------------------------

p3 <- ggplot(
  qc,
  aes(x = percent.mito)
) +
  geom_histogram(
    bins = 100
  ) +
  labs(
    title = "GSE176078 — Mitochondrial RNA Percentage",
    x = "Mitochondrial RNA (%)",
    y = "Number of cells"
  ) +
  theme_classic()

ggsave(
  file.path(qc_dir, "qc_mito_distribution.png"),
  p3,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 11. QC scatter — genes vs counts
# ------------------------------------------------------------

p4 <- ggplot(
  qc,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.4
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Genes Detected vs RNA Counts",
    x = "nCount_RNA",
    y = "nFeature_RNA"
  ) +
  theme_classic()

ggsave(
  file.path(qc_dir, "qc_counts_vs_features.png"),
  p4,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 12. QC scatter — counts vs mitochondrial percentage
# ------------------------------------------------------------

p5 <- ggplot(
  qc,
  aes(
    x = nCount_RNA,
    y = percent.mito
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.4
  ) +
  scale_x_log10() +
  labs(
    title = "RNA Counts vs Mitochondrial Percentage",
    x = "nCount_RNA",
    y = "Mitochondrial RNA (%)"
  ) +
  theme_classic()

ggsave(
  file.path(qc_dir, "qc_counts_vs_mito.png"),
  p5,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 13. QC complete
# ------------------------------------------------------------

message("")
message("==========================================")
message("QC assessment complete")
message("==========================================")
message("Results saved to: ", qc_dir)
message("")
message("No cells were filtered.")
message("QC thresholds will be selected after inspecting")
message("the distributions.")
message("==========================================")