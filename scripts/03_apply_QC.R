# ============================================================
# GSE176078 — Apply QC thresholds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

input_file <- "data/processed/GSE176078_raw_seurat.rds"

qc_dir <- "results/qc"

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Loading Seurat object...")
obj <- readRDS(input_file)

# ------------------------------------------------------------
# QC thresholds
# ------------------------------------------------------------

min_features <- 300
max_features <- 7000

min_counts <- 500
max_counts <- 50000

max_mito <- 15

# ------------------------------------------------------------
# Calculate retention
# ------------------------------------------------------------

meta <- obj@meta.data

keep <- (
  meta$nFeature_RNA > min_features &
  meta$nFeature_RNA < max_features &
  meta$nCount_RNA > min_counts &
  meta$nCount_RNA < max_counts &
  meta$percent.mito < max_mito
)

# ------------------------------------------------------------
# Summary before filtering
# ------------------------------------------------------------

before <- ncol(obj)

message("")
message("==========================================")
message("QC FILTERING")
message("==========================================")
message("Cells before QC: ", before)

message("")
message("Thresholds:")
message("nFeature_RNA: ", min_features, " - ", max_features)
message("nCount_RNA: ", min_counts, " - ", max_counts)
message("percent.mito < ", max_mito)

message("")
message("Cells passing QC: ", sum(keep))
message(
  "Retention: ",
  round(100 * sum(keep) / before, 2),
  "%"
)

message(
  "Removed: ",
  before - sum(keep),
  " cells"
)

# ------------------------------------------------------------
# Breakdown by subtype
# ------------------------------------------------------------

retention_subtype <- data.frame(
  subtype = unique(meta$subtype)
)

retention_subtype$total <- sapply(
  retention_subtype$subtype,
  function(x) sum(meta$subtype == x)
)

retention_subtype$passed <- sapply(
  retention_subtype$subtype,
  function(x) sum(
    meta$subtype == x & keep
  )
)

retention_subtype$retention_percent <- round(
  100 * retention_subtype$passed /
    retention_subtype$total,
  2
)

print(retention_subtype)

fwrite(
  retention_subtype,
  file.path(
    qc_dir,
    "qc_retention_by_subtype.csv"
  )
)

# ------------------------------------------------------------
# Breakdown by patient
# ------------------------------------------------------------

retention_patient <- data.frame(
  patient_id = unique(meta$patient_id)
)

retention_patient$total <- sapply(
  retention_patient$patient_id,
  function(x) sum(meta$patient_id == x)
)

retention_patient$passed <- sapply(
  retention_patient$patient_id,
  function(x) sum(
    meta$patient_id == x & keep
  )
)

retention_patient$retention_percent <- round(
  100 * retention_patient$passed /
    retention_patient$total,
  2
)

retention_patient <- retention_patient[
  order(retention_patient$retention_percent),
]

print(retention_patient)

fwrite(
  retention_patient,
  file.path(
    qc_dir,
    "qc_retention_by_patient.csv"
  )
)

# ------------------------------------------------------------
# Apply filtering
# ------------------------------------------------------------

obj_qc <- subset(
  obj,
  cells = colnames(obj)[keep]
)

# ------------------------------------------------------------
# Save filtered object
# ------------------------------------------------------------

saveRDS(
  obj_qc,
  "data/processed/GSE176078_QC_filtered.rds",
  compress = FALSE
)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

message("")
message("==========================================")
message("QC COMPLETE")
message("==========================================")
message("Final cells: ", ncol(obj_qc))
message("Genes: ", nrow(obj_qc))
message(
  "Patients: ",
  length(unique(obj_qc$patient_id))
)
message("==========================================")