# ============================================================
# GSE113196 — NORMAL BREAST EPITHELIAL PREPROCESSING
# Script 14: Create Seurat object + initial QC
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

message("==============================================")
message("GSE113196 — NORMAL BREAST EPITHELIAL DATA")
message("==============================================")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_dir <- "data/raw/GSE113196/extracted"

output_dir <- "results/normal_breast"

processed_dir <- "data/processed"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Locate matrices
# ------------------------------------------------------------

files <- list.files(
  input_dir,
  pattern = "Expression_Matrix\\.txt\\.gz$",
  full.names = TRUE
)

if (length(files) != 4) {
  stop(
    "Expected 4 GSE113196 expression matrices, found: ",
    length(files)
  )
}

files <- sort(files)

message("")
message("Expression matrices:")
print(basename(files))

# ------------------------------------------------------------
# 3. Function to read one matrix
# ------------------------------------------------------------

read_expression_matrix <- function(file) {

  message("")
  message("Reading: ", basename(file))

  # Read the first column as gene names.
  # data.table::fread can read gzipped text directly.

  dt <- read.delim(
  gzfile(file, open = "rt"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

  message(
    "Raw dimensions: ",
    nrow(dt),
    " rows x ",
    ncol(dt),
    " columns"
  )

  # First column should contain gene identifiers
  gene_column <- colnames(dt)[1]

  genes <- dt[[gene_column]]

  expression_data <- dt[, -1, drop = FALSE]

  # Make sure gene names are unique
  genes <- make.unique(as.character(genes))

  # Convert to numeric matrix
  expression_matrix <- as.matrix(expression_data)

  storage.mode(expression_matrix) <- "numeric"

  rownames(expression_matrix) <- genes

  # Remove completely empty / zero rows
  keep_genes <- Matrix::rowSums(
    expression_matrix != 0
  ) > 0

  expression_matrix <- expression_matrix[keep_genes, , drop = FALSE]

  # Convert to sparse matrix
  expression_matrix <- Matrix(
    expression_matrix,
    sparse = TRUE
  )

  message(
    "Final dimensions: ",
    nrow(expression_matrix),
    " genes x ",
    ncol(expression_matrix),
    " cells"
  )

  return(expression_matrix)
}

# ------------------------------------------------------------
# 4. Read all four individuals
# ------------------------------------------------------------

matrices <- list()

for (file in files) {

  sample_name <- sub(
    "_Expression_Matrix.*",
    "",
    basename(file)
  )

  matrices[[sample_name]] <- read_expression_matrix(file)
}

# ------------------------------------------------------------
# 5. Check dimensions
# ------------------------------------------------------------

message("")
message("==============================================")
message("INDIVIDUAL SAMPLE DIMENSIONS")
message("==============================================")

for (sample_name in names(matrices)) {

  message(
    sample_name,
    ": ",
    nrow(matrices[[sample_name]]),
    " genes x ",
    ncol(matrices[[sample_name]]),
    " cells"
  )
}

# ------------------------------------------------------------
# 6. Standardize gene universe
# ------------------------------------------------------------

message("")
message("Checking gene identifiers...")

common_genes <- Reduce(
  intersect,
  lapply(matrices, rownames)
)

message(
  "Genes shared across all four individuals: ",
  length(common_genes)
)

if (length(common_genes) < 1000) {
  stop(
    "Unexpectedly small number of shared genes. ",
    "Inspect gene identifiers before proceeding."
  )
}

# Use union of genes so no expression information is discarded.
all_genes <- Reduce(
  union,
  lapply(matrices, rownames)
)

message(
  "Total genes across all individuals: ",
  length(all_genes)
)

# ------------------------------------------------------------
# 7. Align matrices to common gene universe
# ------------------------------------------------------------

message("")
message("Aligning gene matrices...")

aligned_matrices <- lapply(
  matrices,
  function(x) {

    missing_genes <- setdiff(
      all_genes,
      rownames(x)
    )

    if (length(missing_genes) > 0) {

      zero_matrix <- Matrix(
        0,
        nrow = length(missing_genes),
        ncol = ncol(x),
        sparse = TRUE
      )

      rownames(zero_matrix) <- missing_genes
      colnames(zero_matrix) <- colnames(x)

      x <- rbind(
        x,
        zero_matrix
      )
    }

    x[all_genes, , drop = FALSE]
  }
)

# ------------------------------------------------------------
# 8. Combine individuals
# ------------------------------------------------------------

message("")
message("Combining individuals...")

combined_counts <- do.call(
  cbind,
  aligned_matrices
)

message(
  "Combined matrix: ",
  nrow(combined_counts),
  " genes x ",
  ncol(combined_counts),
  " cells"
)

# ------------------------------------------------------------
# 9. Ensure unique cell names
# ------------------------------------------------------------

if (anyDuplicated(colnames(combined_counts)) > 0) {

  message("Duplicate cell barcodes detected.")

  colnames(combined_counts) <- make.unique(
    colnames(combined_counts)
  )
}

# ------------------------------------------------------------
# 10. Create sample metadata
# ------------------------------------------------------------

cell_names <- colnames(combined_counts)

sample_id <- sub(
  "_.*",
  "",
  cell_names
)

sample_metadata <- data.frame(
  cell = cell_names,
  sample_id = sample_id,
  individual = sample_id,
  row.names = cell_names,
  stringsAsFactors = FALSE
)

message("")
message("Sample distribution:")
print(table(sample_metadata$sample_id))

# ------------------------------------------------------------
# 11. Create Seurat object
# ------------------------------------------------------------

message("")
message("Creating Seurat object...")

obj <- CreateSeuratObject(
  counts = combined_counts,
  project = "GSE113196_NormalBreast",
  meta.data = sample_metadata,
  min.cells = 3,
  min.features = 100
)

# ------------------------------------------------------------
# 12. QC metrics
# ------------------------------------------------------------

message("")
message("Calculating QC metrics...")

# Human mitochondrial genes
obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# Ribosomal genes
obj[["percent.ribo"]] <- PercentageFeatureSet(
  obj,
  pattern = "^RPS|^RPL"
)

# ------------------------------------------------------------
# 13. QC summary
# ------------------------------------------------------------

message("")
message("==============================================")
message("INITIAL QC SUMMARY")
message("==============================================")

print(
  summary(
    obj@meta.data[
      ,
      c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "percent.ribo"
      )
    ]
  )
)

# ------------------------------------------------------------
# 14. QC by individual
# ------------------------------------------------------------

message("")
message("==============================================")
message("QC BY INDIVIDUAL")
message("==============================================")

qc_by_sample <- obj@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    cells = n(),

    median_genes = median(
      nFeature_RNA
    ),

    median_counts = median(
      nCount_RNA
    ),

    median_mt = median(
      percent.mt
    ),

    median_ribo = median(
      percent.ribo
    ),

    .groups = "drop"
  )

print(qc_by_sample)

write.csv(
  qc_by_sample,
  file.path(
    output_dir,
    "initial_QC_by_individual.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. QC violin plot
# ------------------------------------------------------------

message("")
message("Generating QC violin plot...")

p_qc <- VlnPlot(
  obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.ribo"
  ),
  group.by = "sample_id",
  ncol = 4,
  pt.size = 0
)

ggsave(
  file.path(
    output_dir,
    "initial_QC_violin_by_individual.png"
  ),
  p_qc,
  width = 16,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 16. Save initial object
# ------------------------------------------------------------

output_file <- file.path(
  processed_dir,
  "GSE113196_normal_breast_raw_v1.rds"
)

message("")
message("Saving initial Seurat object...")

saveRDS(
  obj,
  output_file,
  compress = FALSE
)

# ------------------------------------------------------------
# 17. Final summary
# ------------------------------------------------------------

message("")
message("==============================================")
message("SCRIPT 14 COMPLETE")
message("==============================================")

message("")
message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))

message("")
message("Individuals:")
print(table(obj$sample_id))

message("")
message("Saved object:")
message(output_file)

message("")
message("Saved QC report:")
message(
  file.path(
    output_dir,
    "initial_QC_by_individual.csv"
  )
)

message("")
message("Saved QC plot:")
message(
  file.path(
    output_dir,
    "initial_QC_violin_by_individual.png"
  )
)

message("")
message("IMPORTANT:")
message(
  "No cells have been filtered yet. ",
  "QC thresholds will be chosen after inspecting the distributions."
)

message("==============================================")