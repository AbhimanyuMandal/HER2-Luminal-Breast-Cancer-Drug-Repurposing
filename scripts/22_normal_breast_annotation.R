library(Seurat)
library(dplyr)
library(ggplot2)

cat("==============================================\n")
cat("GSE113196 — FINAL NORMAL BREAST ANNOTATION\n")
cat("==============================================\n\n")

# --------------------------------------------------
# Load corrected clustered object
# --------------------------------------------------

x <- readRDS(
  "data/processed/GSE113196_normal_breast_RPCA_UMAP_clustered_v1_fixed.rds"
)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n\n")

# --------------------------------------------------
# Final cluster annotations
# --------------------------------------------------

annotation_map <- c(
  "0"  = "Myoepithelial",
  "1"  = "Luminal_Epithelial",
  "2"  = "Basal_Epithelial",
  "3"  = "Luminal_Epithelial",
  "4"  = "Myoepithelial_Basal",
  "5"  = "Luminal_Secretory",
  "6"  = "Luminal_Secretory",
  "7"  = "Luminal_Epithelial",
  "8"  = "Luminal_Epithelial_Cycling",
  "9"  = "Cycling_Basal_Myoepithelial",
  "10" = "Fibroblast_Stromal",
  "11" = "Unresolved_Singleton"
)

confidence_map <- c(
  "0"  = "High",
  "1"  = "High",
  "2"  = "High",
  "3"  = "High",
  "4"  = "Moderate-High",
  "5"  = "High",
  "6"  = "Moderate-High",
  "7"  = "Moderate",
  "8"  = "Low-Moderate",
  "9"  = "High_Cycling",
  "10" = "Very_High",
  "11" = "Very_Low"
)

# --------------------------------------------------
# Apply annotations
# --------------------------------------------------

cluster_ids <- as.character(x$seurat_clusters)

x$celltype_v1 <- unname(
  annotation_map[cluster_ids]
)

x$annotation_confidence_v1 <- unname(
  confidence_map[cluster_ids]
)

# --------------------------------------------------
# Check
# --------------------------------------------------

cat("\nFINAL ANNOTATION TABLE\n")
print(
  data.frame(
    cluster = names(annotation_map),
    celltype = unname(annotation_map),
    confidence = unname(confidence_map)
  ),
  row.names = FALSE
)

cat("\nCells by cell type:\n")
print(table(x$celltype_v1))

# --------------------------------------------------
# Output directory
# --------------------------------------------------

dir.create(
  "results/normal_breast/annotation",
  recursive = TRUE,
  showWarnings = FALSE
)

# --------------------------------------------------
# Save annotation table
# --------------------------------------------------

annotation_table <- data.frame(
  cluster = names(annotation_map),
  celltype = unname(annotation_map),
  confidence = unname(confidence_map)
)

write.csv(
  annotation_table,
  "results/normal_breast/annotation/final_annotation_map_v1.csv",
  row.names = FALSE
)

# --------------------------------------------------
# Cell counts
# --------------------------------------------------

cell_counts <- x@meta.data %>%
  count(
    celltype_v1,
    name = "cells"
  ) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  )

write.csv(
  cell_counts,
  "results/normal_breast/annotation/final_celltype_counts_v1.csv",
  row.names = FALSE
)

# --------------------------------------------------
# UMAP
# --------------------------------------------------

p1 <- DimPlot(
  x,
  reduction = "umap",
  group.by = "celltype_v1",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("GSE113196 Normal Breast — Final Cell-Type Annotation")

ggsave(
  "results/normal_breast/annotation/UMAP_final_celltypes_v1.png",
  p1,
  width = 12,
  height = 8,
  dpi = 300
)

# --------------------------------------------------
# Confidence UMAP
# --------------------------------------------------

p2 <- DimPlot(
  x,
  reduction = "umap",
  group.by = "annotation_confidence_v1"
) +
  ggtitle("GSE113196 Normal Breast — Annotation Confidence")

ggsave(
  "results/normal_breast/annotation/UMAP_annotation_confidence_v1.png",
  p2,
  width = 12,
  height = 8,
  dpi = 300
)

# --------------------------------------------------
# Save final object
# --------------------------------------------------

saveRDS(
  x,
  "data/processed/GSE113196_normal_breast_annotated_v1.rds"
)

cat("\n==============================================\n")
cat("FINAL ANNOTATION COMPLETE\n")
cat("==============================================\n\n")

cat("Saved:\n")
cat("- data/processed/GSE113196_normal_breast_annotated_v1.rds\n")
cat("- results/normal_breast/annotation/final_annotation_map_v1.csv\n")
cat("- results/normal_breast/annotation/final_celltype_counts_v1.csv\n")
cat("- results/normal_breast/annotation/UMAP_final_celltypes_v1.png\n")
cat("- results/normal_breast/annotation/UMAP_annotation_confidence_v1.png\n")