# ============================================================
# GSE176078 — Subtype Characterization
# Version: V1
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

message("==============================================")
message("SUBTYPE CHARACTERIZATION")
message("==============================================")

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/processed/GSE176078_RPCA_UMAP_annotated_v1.rds"

output_dir <- "results/subtype_characterization"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load annotated object
# ------------------------------------------------------------

message("")
message("Loading annotated object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))
message("Clusters: ", length(unique(Idents(obj))))

# ------------------------------------------------------------
# 3. Check required metadata
# ------------------------------------------------------------

required_metadata <- c(
  "celltype_major_original",
  "celltype_subtype",
  "annotation_confidence"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste(
      "Missing metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 4. Set identity to final subtype
# ------------------------------------------------------------

Idents(obj) <- "celltype_subtype"

# ------------------------------------------------------------
# 5. Subtype abundance
# ------------------------------------------------------------

message("")
message("Calculating subtype abundance...")

subtype_counts <- obj@meta.data %>%
  count(
    celltype_subtype,
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
    "subtype_abundance.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Abundance plot
# ------------------------------------------------------------

p_abundance <- ggplot(
  subtype_counts,
  aes(
    x = reorder(celltype_subtype, cells),
    y = cells
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Cell-Type/Subtype Abundance",
    x = NULL,
    y = "Number of cells"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "subtype_abundance.png"
  ),
  p_abundance,
  width = 10,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Major biological compartments
# ------------------------------------------------------------

message("")
message("Assigning biological compartments...")

compartment_map <- c(

  "Luminal_ERplus_Cancer" = "Tumor_Epithelial",
  "Luminal_Epithelial" = "Tumor_Epithelial",
  "Basal_Epithelial" = "Tumor_Epithelial",
  "Cycling_Epithelial_LowConfidence" = "Tumor_Epithelial",
  "Cycling_Lineage_Unresolved" = "Unresolved_Cycling",

  "NK" = "Immune",
  "CD4_T" = "Immune",
  "Treg" = "Immune",
  "Activated_Cytotoxic_T" = "Immune",
  "B_cell" = "Immune",
  "Plasmablast_Plasma" = "Immune",
  "FOLR2_Resident_Macrophage" = "Immune",
  "FCN1_Inflammatory_Myeloid" = "Immune",
  "Plasmacytoid_Dendritic_Cell" = "Immune",
  "Cycling_Myeloid_Unresolved" = "Unresolved_Cycling",

  "Fibroblast_like_CAF" = "Stromal",
  "LRRC15_Matrix_CAF" = "Stromal",

  "ACKR1_Endothelial" = "Endothelial",
  "Vascular_Endothelial" = "Endothelial",
  "Lymphatic_Endothelial" = "Endothelial",

  "Smooth_Muscle_Perivascular" = "Perivascular",
  "RGS5_Pericyte" = "Perivascular",

  "Cycling_Endothelial_Unresolved" = "Unresolved_Cycling"
)

obj$biological_compartment <- unname(
  compartment_map[as.character(obj$celltype_subtype)]
)

if (any(is.na(obj$biological_compartment))) {
  stop("Some subtypes do not have a biological compartment.")
}

# ------------------------------------------------------------
# 8. Compartment abundance
# ------------------------------------------------------------

compartment_counts <- obj@meta.data %>%
  count(
    biological_compartment,
    name = "cells"
  ) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  ) %>%
  arrange(desc(cells))

write.csv(
  compartment_counts,
  file.path(
    output_dir,
    "biological_compartment_abundance.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Compartment UMAP
# ------------------------------------------------------------

message("")
message("Generating compartment UMAP...")

p_compartment <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = "biological_compartment",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle(
    "Major Biological Compartments"
  )

ggsave(
  file.path(
    output_dir,
    "UMAP_biological_compartments.png"
  ),
  p_compartment,
  width = 11,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Canonical marker sets
# ------------------------------------------------------------

marker_sets <- list(

  # Tumor / epithelial
  Luminal = c(
    "ESR1",
    "PGR",
    "FOXA1",
    "GATA3",
    "TFF1",
    "AGR2",
    "MUC1",
    "KRT8",
    "KRT18",
    "KRT19"
  ),

  Basal = c(
    "KRT5",
    "KRT14",
    "KRT17",
    "KRT15",
    "TP63",
    "COL17A1"
  ),

  # T / NK
  T_cell = c(
    "CD3D",
    "CD3E",
    "CD3G",
    "TRBC1",
    "TRBC2"
  ),

  Cytotoxic = c(
    "NKG7",
    "GNLY",
    "GZMB",
    "GZMH",
    "PRF1",
    "CCL5"
  ),

  Treg = c(
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "TNFRSF4",
    "TIGIT"
  ),

  # B / plasma
  B_cell = c(
    "CD79A",
    "CD79B",
    "MS4A1",
    "CD19",
    "CD74",
    "HLA-DRA"
  ),

  Plasma = c(
    "MZB1",
    "JCHAIN",
    "DERL3",
    "IGHG1",
    "IGKC"
  ),

  # Myeloid
  Macrophage = c(
    "C1QA",
    "C1QB",
    "C1QC",
    "CD68",
    "MS4A7",
    "FOLR2",
    "CD163"
  ),

  Inflammatory_Myeloid = c(
    "FCN1",
    "CD14",
    "S100A8",
    "S100A9",
    "IL1B",
    "TREM1",
    "LYZ"
  ),

  pDC = c(
    "CLEC4C",
    "LILRA4",
    "GZMB",
    "PLD4",
    "SPIB"
  ),

  # CAF
  Fibroblast = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "SPARC",
    "COL6A1"
  ),

  Matrix_CAF = c(
    "LRRC15",
    "COL11A1",
    "COL10A1",
    "POSTN",
    "CTHRC1",
    "MMP11"
  ),

  # Endothelial
  Endothelial = c(
    "PECAM1",
    "VWF",
    "KDR",
    "CDH5",
    "EMCN",
    "CLDN5"
  ),

  ACKR1_Endothelial = c(
    "ACKR1",
    "SELE",
    "SELP",
    "VWF"
  ),

  Lymphatic = c(
    "CCL21",
    "PROX1",
    "FLT4",
    "PDPN",
    "LYVE1"
  ),

  # Perivascular
  Pericyte = c(
    "RGS5",
    "KCNJ8",
    "ABCC9",
    "PDGFRB",
    "COX4I2"
  ),

  Smooth_Muscle = c(
    "ACTA2",
    "MYH11",
    "TAGLN",
    "CNN1",
    "LMOD1"
  ),

  # Cycling
  Cycling = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "BIRC5",
    "CENPF",
    "STMN1",
    "TYMS",
    "RRM2",
    "PLK1"
  )
)

# ------------------------------------------------------------
# 11. DotPlot — all major marker programs
# ------------------------------------------------------------

message("")
message("Generating canonical marker DotPlot...")

all_markers <- unique(
  unlist(marker_sets)
)

available_markers <- intersect(
  all_markers,
  rownames(obj)
)

p_dot <- DotPlot(
  obj,
  features = available_markers,
  group.by = "celltype_subtype",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle(
    "Canonical Marker Expression Across Cell Subtypes"
  )

ggsave(
  file.path(
    output_dir,
    "DotPlot_canonical_markers.png"
  ),
  p_dot,
  width = 18,
  height = 12,
  dpi = 300
)

# ------------------------------------------------------------
# 12. Tumor / epithelial characterization
# ------------------------------------------------------------

message("")
message("Generating epithelial characterization...")

epithelial_subtypes <- c(
  "Luminal_ERplus_Cancer",
  "Luminal_Epithelial",
  "Basal_Epithelial",
  "Cycling_Epithelial_LowConfidence",
  "Cycling_Lineage_Unresolved"
)

epithelial_obj <- subset(
  obj,
  subset = celltype_subtype %in% epithelial_subtypes
)

Idents(epithelial_obj) <- "celltype_subtype"

epithelial_markers <- unique(
  c(
    marker_sets$Luminal,
    marker_sets$Basal,
    marker_sets$Cycling
  )
)

epithelial_markers <- intersect(
  epithelial_markers,
  rownames(epithelial_obj)
)

p_epithelial <- DotPlot(
  epithelial_obj,
  features = epithelial_markers,
  group.by = "celltype_subtype",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle(
    "Epithelial / Tumor Cell Characterization"
  )

ggsave(
  file.path(
    output_dir,
    "DotPlot_epithelial_characterization.png"
  ),
  p_epithelial,
  width = 14,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 13. Immune characterization
# ------------------------------------------------------------

message("")
message("Generating immune characterization...")

immune_subtypes <- c(
  "NK",
  "CD4_T",
  "Treg",
  "Activated_Cytotoxic_T",
  "B_cell",
  "Plasmablast_Plasma",
  "FOLR2_Resident_Macrophage",
  "FCN1_Inflammatory_Myeloid",
  "Plasmacytoid_Dendritic_Cell",
  "Cycling_Myeloid_Unresolved"
)

immune_obj <- subset(
  obj,
  subset = celltype_subtype %in% immune_subtypes
)

Idents(immune_obj) <- "celltype_subtype"

immune_markers <- unique(
  c(
    marker_sets$T_cell,
    marker_sets$Cytotoxic,
    marker_sets$Treg,
    marker_sets$B_cell,
    marker_sets$Plasma,
    marker_sets$Macrophage,
    marker_sets$Inflammatory_Myeloid,
    marker_sets$pDC,
    marker_sets$Cycling
  )
)

immune_markers <- intersect(
  immune_markers,
  rownames(immune_obj)
)

p_immune <- DotPlot(
  immune_obj,
  features = immune_markers,
  group.by = "celltype_subtype",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle(
    "Immune Cell-State Characterization"
  )

ggsave(
  file.path(
    output_dir,
    "DotPlot_immune_characterization.png"
  ),
  p_immune,
  width = 18,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 14. Stromal / vascular characterization
# ------------------------------------------------------------

message("")
message("Generating stromal characterization...")

stromal_subtypes <- c(
  "Fibroblast_like_CAF",
  "LRRC15_Matrix_CAF",
  "ACKR1_Endothelial",
  "Vascular_Endothelial",
  "Lymphatic_Endothelial",
  "Smooth_Muscle_Perivascular",
  "RGS5_Pericyte"
)

stromal_obj <- subset(
  obj,
  subset = celltype_subtype %in% stromal_subtypes
)

Idents(stromal_obj) <- "celltype_subtype"

stromal_markers <- unique(
  c(
    marker_sets$Fibroblast,
    marker_sets$Matrix_CAF,
    marker_sets$Endothelial,
    marker_sets$ACKR1_Endothelial,
    marker_sets$Lymphatic,
    marker_sets$Pericyte,
    marker_sets$Smooth_Muscle
  )
)

stromal_markers <- intersect(
  stromal_markers,
  rownames(stromal_obj)
)

p_stromal <- DotPlot(
  stromal_obj,
  features = stromal_markers,
  group.by = "celltype_subtype",
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle(
    "Stromal and Vascular Characterization"
  )

ggsave(
  file.path(
    output_dir,
    "DotPlot_stromal_characterization.png"
  ),
  p_stromal,
  width = 18,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Feature plots for key biological markers
# ------------------------------------------------------------

message("")
message("Generating key FeaturePlots...")

key_features <- intersect(
  c(
    "ESR1",
    "TFF1",
    "FOXA1",
    "KRT5",
    "KRT14",
    "MKI67",
    "NKG7",
    "FOXP3",
    "FOLR2",
    "LRRC15",
    "RGS5",
    "ACKR1",
    "CCL21"
  ),
  rownames(obj)
)

p_features <- FeaturePlot(
  obj,
  features = key_features,
  reduction = "umap.rpca",
  raster = TRUE,
  ncol = 4
)

ggsave(
  file.path(
    output_dir,
    "FeaturePlot_key_biological_markers.png"
  ),
  p_features,
  width = 16,
  height = 12,
  dpi = 300
)

# ------------------------------------------------------------
# 16. Save updated object with compartment metadata
# ------------------------------------------------------------

output_file <- "data/processed/GSE176078_RPCA_UMAP_annotated_v1_characterized.rds"

message("")
message("Saving characterized object...")

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
message("SUBTYPE CHARACTERIZATION COMPLETE")
message("==============================================")

message("")
message("Biological compartments:")
print(compartment_counts)

message("")
message("Output directory:")
message(output_dir)

message("")
message("Characterized object:")
message(output_file)

message("")
message("==============================================")