# ============================================================
# GSE176078 — Final Marker-Driven Cell-Type Annotation
# Version: V1
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

message("==============================================")
message("FINAL MARKER-DRIVEN ANNOTATION")
message("==============================================")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/processed/GSE176078_RPCA_UMAP_clustered_joined.rds"

output_file <- "data/processed/GSE176078_RPCA_UMAP_annotated_v1.rds"

output_dir <- "results/annotation"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load clustered object
# ------------------------------------------------------------

message("")
message("Loading clustered object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))
message("Clusters: ", length(unique(Idents(obj))))

# ------------------------------------------------------------
# 3. Check required metadata
# ------------------------------------------------------------

if (!"celltype_major" %in% colnames(obj@meta.data)) {
  stop(
    "Required metadata column 'celltype_major' was not found."
  )
}

# Preserve original annotation
obj$celltype_major_original <- obj$celltype_major

# ------------------------------------------------------------
# 4. Final cluster annotations
#
# These labels were determined using:
#   - Cluster marker genes
#   - Marker-program analysis
#   - Original celltype_major composition
#   - Biological marker evidence
#
# Ambiguous cycling clusters are deliberately retained as
# unresolved rather than being over-interpreted.
# ------------------------------------------------------------

final_annotations <- c(

  "0"  = "Luminal_ERplus_Cancer",

  "1"  = "NK",

  "2"  = "CD4_T",

  "3"  = "Treg",

  "4"  = "Activated_Cytotoxic_T",

  "5"  = "FOLR2_Resident_Macrophage",

  "6"  = "Basal_Epithelial",

  "7"  = "Luminal_Epithelial",

  "8"  = "ACKR1_Endothelial",

  "9"  = "Cycling_Lineage_Unresolved",

  "10" = "Plasmablast_Plasma",

  "11" = "Fibroblast_like_CAF",

  "12" = "Vascular_Endothelial",

  "13" = "B_cell",

  "14" = "LRRC15_Matrix_CAF",

  "15" = "Smooth_Muscle_Perivascular",

  "16" = "FCN1_Inflammatory_Myeloid",

  "17" = "RGS5_Pericyte",

  "18" = "Basal_Epithelial",

  "19" = "Cycling_Myeloid_Unresolved",

  "20" = "Plasmacytoid_Dendritic_Cell",

  "21" = "Lymphatic_Endothelial",

  "22" = "Cycling_Epithelial_LowConfidence"
)

# ------------------------------------------------------------
# 5. Annotation confidence
# ------------------------------------------------------------

annotation_confidence <- c(

  "0"  = "High",
  "1"  = "High",
  "2"  = "High",
  "3"  = "High",
  "4"  = "High",
  "5"  = "High",
  "6"  = "High",
  "7"  = "High",
  "8"  = "High",
  "9"  = "Low",
  "10" = "High",
  "11" = "High",
  "12" = "High",
  "13" = "High",
  "14" = "High",
  "15" = "High",
  "16" = "High",
  "17" = "High",
  "18" = "High",
  "19" = "Low",
  "20" = "High",
  "21" = "High",
  "22" = "Low"
)

# ------------------------------------------------------------
# 6. Apply annotations to cells
# ------------------------------------------------------------

cluster_ids <- as.character(Idents(obj))

obj$celltype_subtype <- unname(
  final_annotations[cluster_ids]
)

obj$annotation_confidence <- unname(
  annotation_confidence[cluster_ids]
)

# ------------------------------------------------------------
# 7. Validate annotation mapping
# ------------------------------------------------------------

if (any(is.na(obj$celltype_subtype))) {

  missing_clusters <- unique(
    cluster_ids[is.na(obj$celltype_subtype)]
  )

  stop(
    paste(
      "Missing annotation for cluster(s):",
      paste(missing_clusters, collapse = ", ")
    )
  )
}

if (any(is.na(obj$annotation_confidence))) {

  missing_clusters <- unique(
    cluster_ids[is.na(obj$annotation_confidence)]
  )

  stop(
    paste(
      "Missing confidence for cluster(s):",
      paste(missing_clusters, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 8. Create cluster annotation table
# ------------------------------------------------------------

cluster_annotation_table <- data.frame(
  cluster = names(final_annotations),
  celltype_subtype = unname(final_annotations),
  annotation_confidence = unname(annotation_confidence),
  stringsAsFactors = FALSE
)

# Add original annotation composition

original_annotation <- obj@meta.data %>%
  mutate(
    cluster = as.character(seurat_clusters)
  ) %>%
  count(
    cluster,
    celltype_major_original,
    name = "cells"
  ) %>%
  group_by(cluster) %>%
  mutate(
    total_cluster_cells = sum(cells),
    percentage = 100 * cells / total_cluster_cells
  ) %>%
  arrange(
    cluster,
    desc(percentage)
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    cluster,
    original_major_annotation = celltype_major_original,
    original_annotation_percentage = percentage
  )

cluster_annotation_table <- cluster_annotation_table %>%
  left_join(
    original_annotation,
    by = "cluster"
  ) %>%
  arrange(
    as.numeric(cluster)
  )

# ------------------------------------------------------------
# 9. Save final annotation table
# ------------------------------------------------------------

write.csv(
  cluster_annotation_table,
  file.path(
    output_dir,
    "final_cluster_annotations.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Cell counts by final subtype
# ------------------------------------------------------------

subtype_counts <- obj@meta.data %>%
  count(
    celltype_subtype,
    annotation_confidence,
    name = "cells"
  ) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  ) %>%
  arrange(
    desc(cells)
  )

write.csv(
  subtype_counts,
  file.path(
    output_dir,
    "final_subtype_cell_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Original major annotation × final subtype
# ------------------------------------------------------------

major_subtype_table <- table(
  obj$celltype_major_original,
  obj$celltype_subtype
)

write.csv(
  as.data.frame(major_subtype_table),
  file.path(
    output_dir,
    "original_major_vs_final_subtype.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Cluster sizes
# ------------------------------------------------------------

cluster_sizes <- obj@meta.data %>%
  count(
    seurat_clusters,
    name = "cells"
  ) %>%
  rename(
    cluster = seurat_clusters
  ) %>%
  mutate(
    cluster = as.character(cluster)
  ) %>%
  left_join(
    cluster_annotation_table %>%
      select(
        cluster,
        celltype_subtype,
        annotation_confidence
      ),
    by = "cluster"
  ) %>%
  arrange(
    as.numeric(cluster)
  )

write.csv(
  cluster_sizes,
  file.path(
    output_dir,
    "final_cluster_sizes.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Final annotated UMAP
# ------------------------------------------------------------

message("")
message("Generating final annotated UMAP...")

if (!"umap.rpca" %in% Reductions(obj)) {

  warning(
    "Reduction 'umap.rpca' not found. ",
    "Using default UMAP reduction instead."
  )

  umap_reduction <- "umap"

} else {

  umap_reduction <- "umap.rpca"
}

p_umap <- DimPlot(
  obj,
  reduction = umap_reduction,
  group.by = "celltype_subtype",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle(
    "GSE176078 — Final Cell-Type Annotation"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_final_celltype_annotation.png"
  ),
  plot = p_umap,
  width = 14,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 14. Annotation confidence UMAP
# ------------------------------------------------------------

p_confidence <- DimPlot(
  obj,
  reduction = umap_reduction,
  group.by = "annotation_confidence",
  raster = TRUE
) +
  ggtitle(
    "GSE176078 — Annotation Confidence"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_annotation_confidence.png"
  ),
  plot = p_confidence,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Save final annotated Seurat object
# ------------------------------------------------------------

message("")
message("Saving final annotated object...")

saveRDS(
  obj,
  output_file,
  compress = FALSE
)

# ------------------------------------------------------------
# 16. Print final summary
# ------------------------------------------------------------

message("")
message("==============================================")
message("FINAL ANNOTATION COMPLETE")
message("==============================================")

message("")
message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))
message("Clusters: ", length(unique(Idents(obj))))

message("")
message("FINAL CLUSTER ANNOTATIONS:")
print(
  cluster_annotation_table,
  row.names = FALSE
)

message("")
message("FINAL SUBTYPE COUNTS:")
print(
  subtype_counts)

message("")
message("Files saved:")
message(
  "  ", output_file
)
message(
  "  ",
  file.path(
    output_dir,
    "final_cluster_annotations.csv"
  )
)
message(
  "  ",
  file.path(
    output_dir,
    "final_subtype_cell_counts.csv"
  )
)
message(
  "  ",
  file.path(
    output_dir,
    "original_major_vs_final_subtype.csv"
  )
)
message(
  "  ",
  file.path(
    output_dir,
    "final_cluster_sizes.csv"
  )
)
message(
  "  ",
  file.path(
    output_dir,
    "UMAP_final_celltype_annotation.png"
  )
)
message(
  "  ",
  file.path(
    output_dir,
    "UMAP_annotation_confidence.png"
  )
)

message("")
message("==============================================")
message("V1 ANNOTATION CHECKPOINT CREATED")
message("==============================================")