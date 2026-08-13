# ============================================================
# GSE176078 — Inspect PCA structure
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

input_file <- "data/processed/GSE176078_QC_normalized_PCA.rds"

results_dir <- "results/pca"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

message("Loading PCA object...")

obj <- readRDS(input_file)

# ------------------------------------------------------------
# Top genes for each PC
# ------------------------------------------------------------

message("Finding top PCA-associated genes...")

pc_genes <- list()

for (pc in 1:10) {

  message("PC ", pc)

  pc_genes[[paste0("PC", pc)]] <- head(
    Loadings(obj[["pca"]])[
      order(
        abs(
          Loadings(obj[["pca"]])[, pc]
        ),
        decreasing = TRUE
      ),
      pc
    ],
    20
  )
}

# ------------------------------------------------------------
# Save top genes
# ------------------------------------------------------------

pc_table <- do.call(
  rbind,
  lapply(
    names(pc_genes),
    function(pc_name) {

      genes <- names(pc_genes[[pc_name]])

      data.frame(
        PC = pc_name,
        Rank = seq_along(genes),
        Gene = genes,
        Loading = as.numeric(
          pc_genes[[pc_name]]
        )
      )
    }
  )
)

write.csv(
  pc_table,
  file.path(
    results_dir,
    "PCA_top_genes_PC1_PC10.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

for (pc_name in names(pc_genes)) {

  cat("\n====================================\n")
  cat(pc_name, "\n")
  cat("====================================\n")

  print(
    head(
      pc_genes[[pc_name]],
      15
    )
  )
}

# ------------------------------------------------------------
# PCA plots colored by metadata
# ------------------------------------------------------------

p_subtype <- DimPlot(
  obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "subtype",
  raster = TRUE
) +
  ggtitle("PCA: PC1 vs PC2 — Breast Cancer Subtype")

ggsave(
  file.path(
    results_dir,
    "PCA_PC1_PC2_subtype.png"
  ),
  p_subtype,
  width = 8,
  height = 6,
  dpi = 300
)

p_patient <- DimPlot(
  obj,
  reduction = "pca",
  dims = c(1, 2),
  group.by = "patient_id",
  raster = TRUE
) +
  ggtitle("PCA: PC1 vs PC2 — Patient")

ggsave(
  file.path(
    results_dir,
    "PCA_PC1_PC2_patient.png"
  ),
  p_patient,
  width = 9,
  height = 7,
  dpi = 300
)

message("")
message("PCA inspection complete.")