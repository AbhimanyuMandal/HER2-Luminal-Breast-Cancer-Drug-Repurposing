library(data.table)

message("==============================================")
message("GSE113196 — NORMAL BREAST DATA EXTRACTION")
message("==============================================")

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

tar_file <- "data/raw/GSE113196/GSE113196_RAW.tar"
output_dir <- "data/raw/GSE113196/extracted"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 1. Check TAR archive
# ---------------------------------------------------------

if (!file.exists(tar_file)) {
  stop("GSE113196 TAR file not found: ", tar_file)
}

message("")
message("TAR archive found:")
message(tar_file)

# ---------------------------------------------------------
# 2. List files in archive
# ---------------------------------------------------------

message("")
message("==============================================")
message("FILES IN ARCHIVE")
message("==============================================")

archive_files <- system2(
  "tar",
  args = c("-tf", tar_file),
  stdout = TRUE
)

print(archive_files)

# ---------------------------------------------------------
# 3. Extract archive
# ---------------------------------------------------------

message("")
message("==============================================")
message("EXTRACTING")
message("==============================================")

system2(
  "tar",
  args = c("-xf", tar_file, "-C", output_dir)
)

message("Extraction complete.")

# ---------------------------------------------------------
# 4. Locate expression matrices
# ---------------------------------------------------------

matrix_files <- list.files(
  output_dir,
  pattern = "Expression_Matrix\\.txt\\.gz$",
  full.names = TRUE
)

if (length(matrix_files) == 0) {
  stop("No Expression_Matrix.txt.gz files found after extraction.")
}

message("")
message("==============================================")
message("EXPRESSION MATRICES")
message("==============================================")

print(matrix_files)

# ---------------------------------------------------------
# 5. Inspect each matrix
# ---------------------------------------------------------

message("")
message("==============================================")
message("MATRIX INSPECTION")
message("==============================================")

inspection <- list()

for (file in matrix_files) {

  message("")
  message("----------------------------------------------")
  message("FILE: ", basename(file))
  message("----------------------------------------------")

  con <- gzfile(file, open = "rt")

  lines <- readLines(con, n = 6)

  close(con)

  message("")
  message("First 6 lines:")

  cat(
    paste(lines, collapse = "\n"),
    "\n"
  )

  inspection[[basename(file)]] <- lines
}

# ---------------------------------------------------------
# 6. Save inspection report
# ---------------------------------------------------------

report_file <- file.path(
  "results",
  "normal_breast",
  "GSE113196_matrix_inspection.txt"
)

dir.create(
  dirname(report_file),
  recursive = TRUE,
  showWarnings = FALSE
)

con <- file(report_file, open = "wt")

writeLines(
  "GSE113196 MATRIX INSPECTION",
  con
)

writeLines(
  "==============================================",
  con
)

for (name in names(inspection)) {

  writeLines("", con)
  writeLines("----------------------------------------------", con)
  writeLines(paste("FILE:", name), con)
  writeLines("----------------------------------------------", con)

  writeLines(
    inspection[[name]],
    con
  )
}

close(con)

message("")
message("==============================================")
message("INSPECTION COMPLETE")
message("==============================================")

message("Extracted matrices:")
message(output_dir)

message("")
message("Inspection report:")
message(report_file)

message("")
message("IMPORTANT:")
message("No Seurat object has been created yet.")
message("No QC filtering has been performed.")
message("No normalization has been performed.")
message("No annotation has been assigned.")
message("==============================================")