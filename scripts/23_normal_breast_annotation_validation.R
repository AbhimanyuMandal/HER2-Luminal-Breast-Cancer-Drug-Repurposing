# ============================================================
# 23_normal_breast_annotation_validation.R
# GSE113196 — FINAL NORMAL BREAST ANNOTATION VALIDATION
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

cat("\n")
cat("==============================================\n")
cat("GSE113196 — FINAL NORMAL BREAST VALIDATION\n")
cat("==============================================\n\n")


# ============================================================
# 1. LOAD ANNOTATED OBJECT
# ============================================================

x <- readRDS(
  "data/processed/GSE113196_normal_breast_annotated_v1.rds"
)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n\n")


# ============================================================
# 2. PREPARE RNA ASSAY
# ============================================================

DefaultAssay(x) <- "RNA"

cat("Checking RNA layers...\n")

rna_layers <- Layers(x[["RNA"]])
print(rna_layers)

if (length(rna_layers) > 1) {

  cat("\nJoining RNA layers...\n")

  x[["RNA"]] <- JoinLayers(
    x[["RNA"]]
  )

  cat("RNA layers joined.\n\n")

} else {

  cat("RNA layers already joined.\n\n")

}


# ============================================================
# 3. CHECK REQUIRED METADATA
# ============================================================

if (!"seurat_clusters" %in% colnames(x@meta.data)) {
  stop("ERROR: seurat_clusters column not found.")
}

if (!"celltype_v1" %in% colnames(x@meta.data)) {
  stop("ERROR: celltype_v1 column not found.")
}

if (!"annotation_confidence_v1" %in% colnames(x@meta.data)) {
  stop("ERROR: annotation_confidence_v1 column not found.")
}


# ============================================================
# 4. CANONICAL MARKER PROGRAMS
# ============================================================

marker_panels <- list(

  Luminal_Epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "MUC1"
  ),

  Luminal_Secretory = c(
    "KRT8",
    "KRT18",
    "KRT19",
    "MUC1",
    "SLPI",
    "LCN2",
    "SCGB2A2",
    "KRT7",
    "TFF3",
    "PIP",
    "AGR3",
    "SPDEF"
  ),

  Basal_Epithelial = c(
    "KRT5",
    "KRT14",
    "KRT17",
    "KRT19",
    "KRT6A",
    "KRT6B"
  ),

  Myoepithelial = c(
    "ACTA2",
    "TAGLN",
    "MYL9",
    "TPM2",
    "MYLK",
    "CALD1",
    "CAV1"
  ),

  Fibroblast_Stromal = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "COL6A1",
    "COL6A2",
    "PDGFRA",
    "DPT",
    "SFRP4"
  ),

  Cycling = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "PCNA",
    "TYMS",
    "STMN1",
    "CENPF",
    "AURKB",
    "NDC80",
    "HJURP"
  ),

  Immune = c(
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD74",
    "HLA-DRA",
    "LYZ"
  )
)


# ============================================================
# 5. CHECK MARKER AVAILABILITY
# ============================================================

all_markers <- unique(
  unlist(marker_panels)
)

available_markers <- intersect(
  all_markers,
  rownames(x)
)

cat("Marker genes available:",
    length(available_markers),
    "/",
    length(all_markers),
    "\n\n")

for (program_name in names(marker_panels)) {

  present <- intersect(
    marker_panels[[program_name]],
    rownames(x)
  )

  cat(
    program_name,
    ":",
    length(present),
    "/",
    length(marker_panels[[program_name]]),
    "markers present\n"
  )
}


# ============================================================
# 6. CALCULATE CLUSTER-LEVEL EXPRESSION
# ============================================================

cat("\n==============================================\n")
cat("CALCULATING CLUSTER-LEVEL EXPRESSION\n")
cat("==============================================\n\n")

avg <- AverageExpression(
  x,
  assays = "RNA",
  features = available_markers,
  group.by = "seurat_clusters",
  slot = "data"
)$RNA

avg <- as.matrix(avg)

cat(
  "Average-expression matrix:",
  nrow(avg),
  "genes x",
  ncol(avg),
  "clusters\n\n"
)


# ============================================================
# 7. CALCULATE MARKER PROGRAM SCORES
# ============================================================

program_scores <- data.frame(
  cluster = colnames(avg),
  stringsAsFactors = FALSE
)

for (program_name in names(marker_panels)) {

  genes <- intersect(
    marker_panels[[program_name]],
    rownames(avg)
  )

  if (length(genes) > 0) {

    program_scores[[program_name]] <- colMeans(
      avg[genes, , drop = FALSE],
      na.rm = TRUE
    )

  } else {

    program_scores[[program_name]] <- NA_real_

  }
}


# ------------------------------------------------------------
# Determine strongest program
# ------------------------------------------------------------

program_columns <- names(marker_panels)

program_scores$best_program <- apply(
  program_scores[, program_columns, drop = FALSE],
  1,
  function(z) {

    if (all(is.na(z))) {
      return(NA_character_)
    }

    program_columns[
      which.max(
        replace(z, is.na(z), -Inf)
      )
    ]
  }
)

program_scores$best_score <- apply(
  program_scores[, program_columns, drop = FALSE],
  1,
  function(z) {

    if (all(is.na(z))) {
      return(NA_real_)
    }

    max(
      z,
      na.rm = TRUE
    )
  }
)


# ============================================================
# 8. PRINT PROGRAM SCORES
# ============================================================

cat("\n==============================================\n")
cat("MARKER PROGRAM SCORES\n")
cat("==============================================\n\n")

# Numeric scores only
numeric_cols <- sapply(
  program_scores,
  is.numeric
)

print(
  round(
    program_scores[, numeric_cols, drop = FALSE],
    3
  )
)

cat("\n==============================================\n")
cat("STRONGEST PROGRAM PER CLUSTER\n")
cat("==============================================\n\n")

print(
  program_scores[
    ,
    c("cluster", "best_program", "best_score")
  ],
  row.names = FALSE
)


# ============================================================
# 9. FINAL MANUAL ANNOTATION CORRECTIONS
# ============================================================
#
# These corrections are evidence-based and intentionally
# override the automated/initial annotation where necessary.
#
# Cluster 8:
# Luminal epithelial identity with weaker cycling evidence.
#
# Cluster 9:
# Extremely strong cycling program, therefore kept separate.
#
# Cluster 11:
# Only two cells, therefore unresolved.
# ============================================================

cat("\n==============================================\n")
cat("APPLYING FINAL ANNOTATION CORRECTIONS\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# Cluster 8
# ------------------------------------------------------------

x@meta.data$celltype_v1[
  x@meta.data$seurat_clusters == 8
] <- "Luminal_Epithelial_Cycling"

x@meta.data$annotation_confidence_v1[
  x@meta.data$seurat_clusters == 8
] <- "Low-Moderate"


# ------------------------------------------------------------
# Cluster 9
# ------------------------------------------------------------

x@meta.data$celltype_v1[
  x@meta.data$seurat_clusters == 9
] <- "Cycling_Basal_Myoepithelial"

x@meta.data$annotation_confidence_v1[
  x@meta.data$seurat_clusters == 9
] <- "High_Cycling"


# ------------------------------------------------------------
# Cluster 11
# ------------------------------------------------------------

x@meta.data$celltype_v1[
  x@meta.data$seurat_clusters == 11
] <- "Unresolved_Singleton"

x@meta.data$annotation_confidence_v1[
  x@meta.data$seurat_clusters == 11
] <- "Very_Low"


# ============================================================
# 10. FINAL ANNOTATION MAP
# ============================================================

annotation_map <- x@meta.data %>%
  select(
    seurat_clusters,
    celltype_v1,
    annotation_confidence_v1
  ) %>%
  distinct() %>%
  arrange(
    as.numeric(
      as.character(seurat_clusters)
    )
  )

cat("\n==============================================\n")
cat("FINAL ANNOTATION TABLE\n")
cat("==============================================\n\n")

print(
  annotation_map,
  row.names = FALSE
)


# ============================================================
# 11. FINAL CELL-TYPE COUNTS
# ============================================================

cat("\n==============================================\n")
cat("FINAL CELL-TYPE COUNTS\n")
cat("==============================================\n\n")

celltype_counts <- as.data.frame(
  table(
    x@meta.data$celltype_v1
  )
)

colnames(celltype_counts) <- c(
  "celltype",
  "cells"
)

celltype_counts$percentage <- round(
  100 *
    celltype_counts$cells /
    sum(celltype_counts$cells),
  2
)

print(
  celltype_counts,
  row.names = FALSE
)


# ============================================================
# 12. OUTPUT DIRECTORY
# ============================================================

outdir <- "results/normal_breast/annotation_validation"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 13. SAVE PROGRAM SCORES
# ============================================================

write.csv(
  program_scores,
  file.path(
    outdir,
    "final_program_scores_all_clusters.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SAVE ANNOTATION MAP
# ============================================================

write.csv(
  annotation_map,
  file.path(
    outdir,
    "final_annotation_map_v1.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. SAVE CELL-TYPE COUNTS
# ============================================================

write.csv(
  celltype_counts,
  file.path(
    outdir,
    "final_celltype_counts_v1.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. FINAL CELL-TYPE UMAP
# ============================================================

cat("\nGenerating final cell-type UMAP...\n")

p1 <- DimPlot(
  x,
  reduction = "umap",
  group.by = "celltype_v1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "GSE113196 Normal Breast — Final Cell-Type Annotation"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  file.path(
    outdir,
    "UMAP_final_celltypes_v1.png"
  ),
  p1,
  width = 14,
  height = 10,
  dpi = 300
)


# ============================================================
# 17. CONFIDENCE UMAP
# ============================================================

cat("Generating annotation-confidence UMAP...\n")

p2 <- DimPlot(
  x,
  reduction = "umap",
  group.by = "annotation_confidence_v1",
  label = FALSE
) +
  ggtitle(
    "GSE113196 Normal Breast — Annotation Confidence"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  file.path(
    outdir,
    "UMAP_annotation_confidence_v1.png"
  ),
  p2,
  width = 14,
  height = 10,
  dpi = 300
)


# ============================================================
# 18. SAVE FINAL OBJECT
# ============================================================

saveRDS(
  x,
  "data/processed/GSE113196_normal_breast_final_v1.rds"
)


# ============================================================
# 19. FINAL SUMMARY
# ============================================================

cat("\n")
cat("==============================================\n")
cat("SCRIPT 23 COMPLETE\n")
cat("==============================================\n\n")

cat("Final object:\n")
cat(
  "data/processed/GSE113196_normal_breast_final_v1.rds\n\n"
)

cat("Annotation map:\n")
cat(
  "results/normal_breast/annotation_validation/",
  "final_annotation_map_v1.csv\n\n",
  sep = ""
)

cat("Cell-type counts:\n")
cat(
  "results/normal_breast/annotation_validation/",
  "final_celltype_counts_v1.csv\n\n",
  sep = ""
)

cat("UMAP:\n")
cat(
  "results/normal_breast/annotation_validation/",
  "UMAP_final_celltypes_v1.png\n\n",
  sep = ""
)

cat("Confidence UMAP:\n")
cat(
  "results/normal_breast/annotation_validation/",
  "UMAP_annotation_confidence_v1.png\n"
)

cat("\n==============================================\n")
cat("NORMAL BREAST REFERENCE READY\n")
cat("==============================================\n")