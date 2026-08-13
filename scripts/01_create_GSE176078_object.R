# ============================================================
# GSE176078 — Create Seurat Object
# 100,064 cells | 26 breast cancer patients
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

set.seed(12345)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

base_dir <- "data/raw/GSE176078/Wu_etal_2021_BRCA_scRNASeq"

mtx_file      <- file.path(base_dir, "count_matrix_sparse.mtx")
genes_file    <- file.path(base_dir, "count_matrix_genes.tsv")
barcodes_file <- file.path(base_dir, "count_matrix_barcodes.tsv")
metadata_file <- file.path(base_dir, "metadata.csv")

# ------------------------------------------------------------
# 2. Check that files exist
# ------------------------------------------------------------

required_files <- c(
  mtx_file,
  genes_file,
  barcodes_file,
  metadata_file
)

if (!all(file.exists(required_files))) {
  missing_files <- required_files[!file.exists(required_files)]

  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("All input files found.")

# ------------------------------------------------------------
# 3. Load metadata
# ------------------------------------------------------------

message("Loading metadata...")

metadata <- fread(
  metadata_file,
  data.table = FALSE
)

message(
  "Metadata dimensions: ",
  nrow(metadata),
  " cells × ",
  ncol(metadata),
  " columns"
)

# First column is the cell barcode
cell_barcodes <- as.character(metadata[[1]])

rownames(metadata) <- cell_barcodes

# ------------------------------------------------------------
# 4. Load sparse expression matrix
# ------------------------------------------------------------

message("Loading sparse expression matrix...")
message("This may take several minutes on an 8-GB Mac.")

counts <- ReadMtx(
  mtx = mtx_file,
  features = genes_file,
  cells = barcodes_file,
  feature.column = 1,
  cell.column = 1,
  unique.features = TRUE
)

message(
  "Expression matrix: ",
  nrow(counts),
  " genes × ",
  ncol(counts),
  " cells"
)

# ------------------------------------------------------------
# 5. Validate matrix ↔ metadata
# ------------------------------------------------------------

message("Validating cell barcodes...")

matrix_barcodes <- colnames(counts)
metadata_barcodes <- rownames(metadata)

if (length(matrix_barcodes) != length(metadata_barcodes)) {
  stop("Number of matrix cells and metadata rows do not match.")
}

if (!identical(matrix_barcodes, metadata_barcodes)) {

  message("Barcode order differs — aligning metadata to matrix.")

  if (!setequal(matrix_barcodes, metadata_barcodes)) {
    stop(
      "Matrix and metadata contain different cell barcodes."
    )
  }

  metadata <- metadata[matrix_barcodes, , drop = FALSE]
}

message("Barcode validation passed.")

# ------------------------------------------------------------
# 6. Create Seurat object
# ------------------------------------------------------------

message("Creating Seurat object...")

obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE176078",
  min.cells = 3,
  min.features = 100,
  meta.data = metadata
)

# ------------------------------------------------------------
# 7. Preserve original metadata
# ------------------------------------------------------------

obj$patient_id <- obj$orig.ident

# ------------------------------------------------------------
# 8. Calculate mitochondrial percentage independently
# ------------------------------------------------------------

message("Calculating mitochondrial percentage...")

obj[["percent.mt.calc"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# ------------------------------------------------------------
# 9. Basic validation
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE176078 Seurat object")
message("============================================")

message("Genes:   ", nrow(obj))
message("Cells:   ", ncol(obj))
message("Patients: ", length(unique(obj$patient_id)))

message("")
message("Clinical subtype:")
print(table(obj$subtype))

message("")
message("Major cell types:")
print(table(obj$celltype_major))

message("============================================")

# ------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

message("Saving Seurat object...")

saveRDS(
  obj,
  "data/processed/GSE176078_raw_seurat.rds",
  compress = FALSE
)

message("")
message("DONE.")
message(
  "Saved: data/processed/GSE176078_raw_seurat.rds"
)