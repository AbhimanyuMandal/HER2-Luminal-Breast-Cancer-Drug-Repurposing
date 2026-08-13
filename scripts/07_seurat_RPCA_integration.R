# ============================================================
# GSE176078 — Seurat RPCA Integration
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(ggplot2)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

options(future.globals.maxSize = 6 * 1024^3)

# Important for an 8-GB Mac:
# avoid parallel workers duplicating large Seurat objects
plan("sequential")

input_file <- "data/processed/GSE176078_QC_normalized_PCA.rds"

output_dir <- "data/processed"
results_dir <- "results/integration"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

message("Loading QC/PCA object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))

# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------

if (!"patient_id" %in% colnames(obj[[]])) {
  stop("patient_id is missing from metadata.")
}

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay not found.")
}

if (!"pca" %in% Reductions(obj)) {
  stop("PCA reduction not found.")
}

message("")
message("Patients: ", length(unique(obj$patient_id)))
message("Subtypes: ", paste(unique(obj$subtype), collapse = ", "))
message("")

# ------------------------------------------------------------
# Check existing RNA layers
# ------------------------------------------------------------

message("RNA layers before splitting:")

print(Layers(obj[["RNA"]]))

# ------------------------------------------------------------
# Split RNA assay by patient
# ------------------------------------------------------------

message("")
message("Splitting RNA assay by patient...")

obj[["RNA"]] <- split(
  obj[["RNA"]],
  f = obj$patient_id
)

message("RNA layers after splitting:")

print(Layers(obj[["RNA"]]))

# ------------------------------------------------------------
# Seurat RPCA integration
# ------------------------------------------------------------

message("")
message("==============================================")
message("Starting Seurat RPCA integration")
message("==============================================")
message("")

obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  k.weight = 50,
  verbose = TRUE
)

# ------------------------------------------------------------
# Check integration
# ------------------------------------------------------------

if (!"integrated.rpca" %in% Reductions(obj)) {
  stop("integrated.rpca reduction was not created.")
}

message("")
message("==============================================")
message("RPCA INTEGRATION COMPLETE")
message("==============================================")

message("Reductions:")
print(Reductions(obj))

# ------------------------------------------------------------
# Save integrated object
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "GSE176078_RPCA_integrated.rds"
)

message("")
message("Saving integrated object...")

saveRDS(
  obj,
  output_file,
  compress = FALSE
)

message("")
message("Saved:")
message(output_file)