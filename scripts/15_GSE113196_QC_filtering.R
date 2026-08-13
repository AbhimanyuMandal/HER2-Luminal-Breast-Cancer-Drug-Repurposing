# ============================================================
# GSE113196 — NORMAL BREAST EPITHELIAL DATA
# Script 15: QC Filtering
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

message("==============================================")
message("GSE113196 — QC FILTERING")
message("==============================================")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/processed/GSE113196_normal_breast_raw_v1.rds"

output_dir <- "results/normal_breast"

processed_dir <- "data/processed"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load object
# ------------------------------------------------------------

message("")
message("Loading raw normal breast object...")

obj <- readRDS(input_file)

message("Cells before QC: ", ncol(obj))
message("Genes: ", nrow(obj))

# ------------------------------------------------------------
# 3. QC thresholds
# ------------------------------------------------------------

min_features <- 500
max_features <- 6000
min_counts <- 2000

message("")
message("QC thresholds:")
message("nFeature_RNA >= ", min_features)
message("nFeature_RNA <= ", max_features)
message("nCount_RNA >= ", min_counts)

# ------------------------------------------------------------
# 4. Calculate QC flags
# ------------------------------------------------------------

obj$qc_pass <- with(
  obj@meta.data,
  nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    nCount_RNA >= min_counts
)

# ------------------------------------------------------------
# 5. QC summary
# ------------------------------------------------------------

message("")
message("==============================================")
message("QC SUMMARY")
message("==============================================")

total_cells <- ncol(obj)

retained_cells <- sum(
  obj$qc_pass
)

removed_cells <- total_cells - retained_cells

message("Total cells: ", total_cells)
message("Retained cells: ", retained_cells)
message("Removed cells: ", removed_cells)
message(
  "Overall retention: ",
  round(100 * retained_cells / total_cells, 2),
  "%"
)

# ------------------------------------------------------------
# 6. Removal reasons
# ------------------------------------------------------------

qc_reason <- rep(
  "PASS",
  total_cells
)

md <- obj@meta.data

qc_reason[
  md$nFeature_RNA < min_features
] <- "Low_features"

qc_reason[
  md$nFeature_RNA > max_features
] <- "High_features"

qc_reason[
  md$nCount_RNA < min_counts
] <- "Low_counts"

# If a cell satisfies multiple conditions,
# explicitly identify it.

multiple_fail <- (
  rowSums(
    cbind(
      md$nFeature_RNA < min_features,
      md$nFeature_RNA > max_features,
      md$nCount_RNA < min_counts
    )
  ) > 1
)

qc_reason[multiple_fail] <- "Multiple_QC_failures"

qc_summary <- data.frame(
  qc_category = names(
    table(qc_reason)
  ),
  cells = as.integer(
    table(qc_reason)
  )
)

qc_summary$percentage <- (
  100 *
    qc_summary$cells /
    total_cells
)

print(qc_summary)

write.csv(
  qc_summary,
  file.path(
    output_dir,
    "QC_filtering_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. QC summary by individual
# ------------------------------------------------------------

message("")
message("==============================================")
message("QC BY INDIVIDUAL")
message("==============================================")

qc_by_sample <- obj@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    cells_before = n(),
    cells_retained = sum(qc_pass),
    cells_removed = sum(!qc_pass),
    retention_percent =
      100 * cells_retained / cells_before,
    median_genes =
      median(nFeature_RNA),
    median_counts =
      median(nCount_RNA),
    .groups = "drop"
  )

print(qc_by_sample)

write.csv(
  qc_by_sample,
  file.path(
    output_dir,
    "QC_filtering_by_individual.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. QC distributions before filtering
# ------------------------------------------------------------

message("")
message("Generating pre/post QC visualization...")

md_plot <- obj@meta.data %>%
  select(
    sample_id,
    nFeature_RNA,
    nCount_RNA,
    qc_pass
  )

p_features <- ggplot(
  md_plot,
  aes(
    x = sample_id,
    y = nFeature_RNA
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_hline(
    yintercept = min_features,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = max_features,
    linetype = "dashed"
  ) +
  labs(
    title = "Genes Detected per Cell",
    x = "Individual",
    y = "nFeature_RNA"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "QC_nFeature_by_individual.png"
  ),
  p_features,
  width = 9,
  height = 6,
  dpi = 300
)

p_counts <- ggplot(
  md_plot,
  aes(
    x = sample_id,
    y = nCount_RNA
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_hline(
    yintercept = min_counts,
    linetype = "dashed"
  ) +
  labs(
    title = "UMI Counts per Cell",
    x = "Individual",
    y = "nCount_RNA"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "QC_nCount_by_individual.png"
  ),
  p_counts,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Filter cells
# ------------------------------------------------------------

message("")
message("Applying QC filters...")

obj_filtered <- subset(
  obj,
  subset =
    nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    nCount_RNA >= min_counts
)

# ------------------------------------------------------------
# 10. Remove temporary QC flag from final object
# ------------------------------------------------------------

obj_filtered$qc_pass <- NULL

# ------------------------------------------------------------
# 11. Final QC summary
# ------------------------------------------------------------

message("")
message("==============================================")
message("POST-QC SUMMARY")
message("==============================================")

message(
  "Cells before QC: ",
  ncol(obj)
)

message(
  "Cells after QC: ",
  ncol(obj_filtered)
)

message(
  "Cells removed: ",
  ncol(obj) - ncol(obj_filtered)
)

message(
  "Retention: ",
  round(
    100 * ncol(obj_filtered) / ncol(obj),
    2
  ),
  "%"
)

message("")
message("Cells by individual after QC:")
print(
  table(
    obj_filtered$sample_id
  )
)

# ------------------------------------------------------------
# 12. Save filtered object
# ------------------------------------------------------------

output_file <- file.path(
  processed_dir,
  "GSE113196_normal_breast_QC_v1.rds"
)

message("")
message("Saving filtered object...")

saveRDS(
  obj_filtered,
  output_file,
  compress = FALSE
)

# ------------------------------------------------------------
# 13. Final message
# ------------------------------------------------------------

message("")
message("==============================================")
message("SCRIPT 15 COMPLETE")
message("==============================================")

message("")
message("Saved object:")
message(output_file)

message("")
message("QC reports:")
message(
  file.path(
    output_dir,
    "QC_filtering_summary.csv"
  )
)

message(
  file.path(
    output_dir,
    "QC_filtering_by_individual.csv"
  )
)

message("")
message("QC plots:")
message(
  file.path(
    output_dir,
    "QC_nFeature_by_individual.png"
  )
)

message(
  file.path(
    output_dir,
    "QC_nCount_by_individual.png"
  )
)

message("")
message("==============================================")