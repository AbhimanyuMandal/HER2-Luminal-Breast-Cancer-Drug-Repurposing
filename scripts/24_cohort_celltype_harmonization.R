# 24_cohort_celltype_harmonization.R
#
# GSE176078 BREAST CANCER SUBTYPES vs NORMAL BREAST
# Cohort + cell-type harmonization
#
# Tumor cohort:
#   GSE176078
#   ER+
#   HER2+
#   TNBC
#
# Normal reference:
#   GSE113196
#
# Input:
#   data/processed/GSE176078_RPCA_UMAP_annotated_v1_characterized.rds
#   data/processed/GSE113196_normal_breast_final_v1.rds
#
# Output:
#   results/cohort_harmonization/
#
# Purpose:
#   1. Preserve the original GSE176078 clinical subtype
#   2. Harmonize tumor and normal cell-type annotations
#   3. Compare cell-type composition overall
#   4. Compare ER+, HER2+, and TNBC separately against normal
#   5. Create comparison-ready metadata and objects
#
# IMPORTANT:
#   - This is NOT a re-annotation script.
#   - This is NOT a differential-expression script.
#   - Clinical subtype is patient-level metadata inherited from GSE176078.
#   - Cycling cells are retained as a state and are not forced into a
#     specific lineage.
#   - Script 25/26 should perform downstream disease/subtype analyses.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat("==============================================\n")
cat("SCRIPT 24 — COHORT + CELL-TYPE HARMONIZATION\n")
cat("==============================================\n\n")


# ============================================================
# 1. FILE PATHS
# ============================================================

tumor_file <- "data/processed/GSE176078_RPCA_UMAP_annotated_v1_characterized.rds"

normal_file <- "data/processed/GSE113196_normal_breast_final_v1.rds"

output_dir <- "results/cohort_harmonization"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. CHECK INPUT FILES
# ============================================================

if (!file.exists(tumor_file)) {
  stop(
    "GSE176078 tumor object not found: ",
    tumor_file
  )
}

if (!file.exists(normal_file)) {
  stop(
    "GSE113196 normal breast object not found: ",
    normal_file
  )
}


# ============================================================
# 3. LOAD OBJECTS
# ============================================================

cat("Loading GSE176078 breast cancer object...\n")

tumor <- readRDS(tumor_file)

cat(
  "Tumor cells:",
  ncol(tumor),
  "\n"
)

cat(
  "Tumor genes:",
  nrow(tumor),
  "\n\n"
)


cat("Loading GSE113196 normal breast object...\n")

normal <- readRDS(normal_file)

cat(
  "Normal cells:",
  ncol(normal),
  "\n"
)

cat(
  "Normal genes:",
  nrow(normal),
  "\n\n"
)


# ============================================================
# 4. CHECK REQUIRED METADATA
# ============================================================

required_tumor_metadata <- c(
  "patient_id",
  "subtype",
  "celltype_subtype"
)

missing_tumor <- setdiff(
  required_tumor_metadata,
  colnames(tumor@meta.data)
)

if (length(missing_tumor) > 0) {
  stop(
    "Missing GSE176078 metadata: ",
    paste(missing_tumor, collapse = ", ")
  )
}


required_normal_metadata <- c(
  "celltype_v1"
)

missing_normal <- setdiff(
  required_normal_metadata,
  colnames(normal@meta.data)
)

if (length(missing_normal) > 0) {
  stop(
    "Missing GSE113196 metadata: ",
    paste(missing_normal, collapse = ", ")
  )
}


# ============================================================
# 5. PRESERVE ORIGINAL CLINICAL SUBTYPE
# ============================================================

tumor$clinical_subtype_original <- as.character(
  tumor$subtype
)

cat("==============================================\n")
cat("ORIGINAL GSE176078 CLINICAL SUBTYPE\n")
cat("==============================================\n")

print(
  table(
    tumor$clinical_subtype_original,
    useNA = "ifany"
  )
)

cat("\n")


# ============================================================
# 6. STANDARDIZE CLINICAL SUBTYPE LABELS
# ============================================================
#
# The original subtype values are retained unchanged.
# clinical_subtype is a normalized analysis label.
#
# We deliberately STOP on unknown labels instead of silently
# guessing.
# ============================================================

normalize_clinical_subtype <- function(x) {

  # GSE176078 contains these clinical subtype labels:
  # ER+, HER2+, TNBC
  #
  # Preserve the original labels and map them explicitly.
  # Do NOT remove the "+" symbol.

  x <- trimws(as.character(x))

  result <- dplyr::case_when(
    x == "ER+"   ~ "ER_Positive",
    x == "HER2+" ~ "HER2_Positive",
    x == "TNBC"  ~ "TNBC",
    TRUE         ~ NA_character_
  )

  result
}

tumor$clinical_subtype <- normalize_clinical_subtype(
  tumor$clinical_subtype_original
)


unknown_subtypes <- unique(
  tumor$clinical_subtype_original[
    is.na(tumor$clinical_subtype)
  ]
)

if (length(unknown_subtypes) > 0) {

  stop(
    paste(
      "Unrecognized GSE176078 clinical subtype label(s):",
      paste(unknown_subtypes, collapse = ", ")
    )
  )
}


cat("NORMALIZED CLINICAL SUBTYPE\n")
print(
  table(
    tumor$clinical_subtype,
    useNA = "ifany"
  )
)

cat("\n")


# ============================================================
# 7. VALIDATE PATIENT-LEVEL SUBTYPE CONSISTENCY
# ============================================================

patient_subtype <- tumor@meta.data %>%
  select(
    patient_id,
    clinical_subtype_original,
    clinical_subtype
  ) %>%
  distinct()


patient_subtype_check <- patient_subtype %>%
  count(
    patient_id,
    name = "n_subtypes"
  )


if (any(patient_subtype_check$n_subtypes > 1)) {

  inconsistent_patients <- patient_subtype_check$patient_id[
    patient_subtype_check$n_subtypes > 1
  ]

  stop(
    paste(
      "A patient has multiple clinical subtype labels:",
      paste(inconsistent_patients, collapse = ", ")
    )
  )
}


write.csv(
  patient_subtype,
  file.path(
    output_dir,
    "GSE176078_patient_clinical_subtypes.csv"
  ),
  row.names = FALSE
)


cat("==============================================\n")
cat("NORMALIZED CLINICAL SUBTYPE\n")
cat("==============================================\n")

print(
  table(
    tumor$clinical_subtype
  )
)

cat("\n")


# ============================================================
# 8. CREATE SHARED CELL-TYPE FRAMEWORK
# ============================================================
#
# This is a comparison framework.
#
# We do NOT overwrite the original annotations.
# We deliberately use broad biological categories so that
# tumor and normal datasets can be compared without pretending
# that their annotation schemes are identical.
# ============================================================

standardize_celltype <- function(x) {

  x <- as.character(x)

  result <- rep(
    "Other_Unresolved",
    length(x)
  )

  # ----------------------------------------------------------
  # Cycling
  # ----------------------------------------------------------

  result[
    grepl(
      "cycling",
      x,
      ignore.case = TRUE
    )
  ] <- "Cycling"


  # ----------------------------------------------------------
  # Immune
  # ----------------------------------------------------------

  result[
    grepl(
      paste(
        c(
          "immune",
          "T.cell",
          "T_cell",
          "B.cell",
          "B_cell",
          "macrophage",
          "monocyte",
          "myeloid",
          "dendritic",
          "NK",
          "lymph",
          "plasmablast",
          "plasma",
          "Treg"
        ),
        collapse = "|"
      ),
      x,
      ignore.case = TRUE
    )
  ] <- "Immune"


  # ----------------------------------------------------------
  # Fibroblast / stromal / CAF
  # ----------------------------------------------------------

  result[
    grepl(
      "fibro|stromal|CAF",
      x,
      ignore.case = TRUE
    )
  ] <- "Fibroblast_Stromal"


  # ----------------------------------------------------------
  # Endothelial
  # ----------------------------------------------------------

  result[
    grepl(
      "endothelial",
      x,
      ignore.case = TRUE
    )
  ] <- "Endothelial"


  # ----------------------------------------------------------
  # Perivascular / pericyte / smooth muscle
  # ----------------------------------------------------------

  result[
    grepl(
      "pericyte|perivascular|smooth.?muscle",
      x,
      ignore.case = TRUE
    )
  ] <- "Perivascular"


  # ----------------------------------------------------------
  # Luminal secretory
  # ----------------------------------------------------------

  result[
    grepl(
      "luminal.*secretory|secretory",
      x,
      ignore.case = TRUE
    )
  ] <- "Luminal_Secretory"


  # ----------------------------------------------------------
  # Luminal epithelial
  # ----------------------------------------------------------

  result[
    grepl(
      "luminal.*epithelial|luminal",
      x,
      ignore.case = TRUE
    )
  ] <- "Luminal_Epithelial"


  # ----------------------------------------------------------
  # Basal epithelial
  # ----------------------------------------------------------

  result[
    grepl(
      "basal",
      x,
      ignore.case = TRUE
    )
  ] <- "Basal_Epithelial"


  # ----------------------------------------------------------
  # Myoepithelial
  # ----------------------------------------------------------

  result[
    grepl(
      "myoepithelial",
      x,
      ignore.case = TRUE
    )
  ] <- "Myoepithelial"


  result
}


# ============================================================
# 9. APPLY SHARED CELL-TYPE LABELS
# ============================================================

tumor$original_celltype <- as.character(
  tumor$celltype_subtype
)

normal$original_celltype <- as.character(
  normal$celltype_v1
)


tumor$comparison_celltype <- standardize_celltype(
  tumor$original_celltype
)

normal$comparison_celltype <- standardize_celltype(
  normal$original_celltype
)


# ============================================================
# 10. PRINT STANDARDIZED CELL TYPES
# ============================================================

cat("==============================================\n")
cat("STANDARDIZED GSE176078 CELL TYPES\n")
cat("==============================================\n")

print(
  table(
    tumor$comparison_celltype,
    useNA = "ifany"
  )
)

cat("\n==============================================\n")
cat("STANDARDIZED GSE113196 CELL TYPES\n")
cat("==============================================\n")

print(
  table(
    normal$comparison_celltype,
    useNA = "ifany"
  )
)

cat("\n")


# ============================================================
# 11. FLAG UNRESOLVED CELL TYPES
# ============================================================
#
# We report these rather than silently discarding them.
# ============================================================

unresolved_tumor <- tumor@meta.data %>%
  count(
    original_celltype,
    comparison_celltype,
    name = "cells"
  ) %>%
  filter(
    comparison_celltype == "Other_Unresolved"
  )

unresolved_normal <- normal@meta.data %>%
  count(
    original_celltype,
    comparison_celltype,
    name = "cells"
  ) %>%
  filter(
    comparison_celltype == "Other_Unresolved"
  )


write.csv(
  unresolved_tumor,
  file.path(
    output_dir,
    "unresolved_tumor_celltypes.csv"
  ),
  row.names = FALSE
)

write.csv(
  unresolved_normal,
  file.path(
    output_dir,
    "unresolved_normal_celltypes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. ADD DATASET-LEVEL METADATA
# ============================================================

tumor$dataset <- "Breast_Cancer"

normal$dataset <- "Normal_Breast"

normal$clinical_subtype_original <- "Normal"

normal$clinical_subtype <- factor(
  "Normal",
  levels = c(
    "Normal",
    "ER_Positive",
    "HER2_Positive",
    "TNBC"
  )
)


tumor$clinical_subtype <- factor(
  as.character(tumor$clinical_subtype),
  levels = c(
    "Normal",
    "ER_Positive",
    "HER2_Positive",
    "TNBC"
  )
)


# ============================================================
# 13. CREATE COMPARISON METADATA
# ============================================================

tumor_metadata <- tumor@meta.data %>%
  mutate(
    dataset = "Breast_Cancer",
    clinical_group = as.character(
      clinical_subtype
    ),
    patient_or_sample_id = as.character(
      patient_id
    )
  ) %>%
  select(
    dataset,
    clinical_group,
    patient_or_sample_id,
    original_celltype,
    comparison_celltype,
    everything()
  )


normal_metadata <- normal@meta.data %>%
  mutate(
    dataset = "Normal_Breast",
    clinical_group = "Normal",
    patient_or_sample_id = ifelse(
      "sample_id" %in% colnames(normal@meta.data),
      as.character(sample_id),
      NA_character_
    )
  ) %>%
  select(
    dataset,
    clinical_group,
    patient_or_sample_id,
    original_celltype,
    comparison_celltype,
    everything()
  )


combined_metadata <- bind_rows(
  tumor_metadata,
  normal_metadata
)


write.csv(
  combined_metadata,
  file.path(
    output_dir,
    "combined_cell_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. CELL-TYPE COMPOSITION
# ============================================================

composition_counts <- combined_metadata %>%
  count(
    clinical_group,
    comparison_celltype,
    name = "cells"
  ) %>%
  group_by(
    clinical_group
  ) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  ) %>%
  ungroup()


write.csv(
  composition_counts,
  file.path(
    output_dir,
    "celltype_composition_all_groups.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. WIDE COMPOSITION TABLE
# ============================================================

composition_wide <- composition_counts %>%
  select(
    comparison_celltype,
    clinical_group,
    percentage
  ) %>%
  pivot_wider(
    names_from = clinical_group,
    values_from = percentage,
    values_fill = 0
  )


write.csv(
  composition_wide,
  file.path(
    output_dir,
    "celltype_composition_wide.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SUBTYPE-SPECIFIC DIFFERENCES FROM NORMAL
# ============================================================

normal_reference <- composition_counts %>%
  filter(
    clinical_group == "Normal"
  ) %>%
  select(
    comparison_celltype,
    normal_percentage = percentage
  )


subtype_differences <- composition_counts %>%
  filter(
    clinical_group %in% c(
      "ER_Positive",
      "HER2_Positive",
      "TNBC"
    )
  ) %>%
  left_join(
    normal_reference,
    by = "comparison_celltype"
  ) %>%
  mutate(
    normal_percentage = replace_na(
      normal_percentage,
      0
    ),
    percentage_difference =
      percentage - normal_percentage
  )


write.csv(
  subtype_differences,
  file.path(
    output_dir,
    "celltype_composition_subtype_vs_normal.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. OVERALL COMPOSITION PLOT
# ============================================================

cat("Generating overall composition plot...\n")

p1 <- ggplot(
  composition_counts,
  aes(
    x = clinical_group,
    y = percentage,
    fill = comparison_celltype
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(x * 100),
        "%"
      )
    }
  ) +
  labs(
    title =
      "Cell-Type Composition Across Breast Cancer Subtypes and Normal Breast",
    x = NULL,
    y = "Cell proportion",
    fill = "Comparison cell type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "celltype_composition_all_groups.png"
  ),
  p1,
  width = 11,
  height = 7,
  dpi = 300
)


# ============================================================
# 18. SUBTYPE VS NORMAL COMPOSITION PLOT
# ============================================================

cat("Generating subtype-vs-normal composition plot...\n")

plot_data <- composition_counts %>%
  filter(
    clinical_group %in% c(
      "Normal",
      "ER_Positive",
      "HER2_Positive",
      "TNBC"
    )
  )


p2 <- ggplot(
  plot_data,
  aes(
    x = clinical_group,
    y = percentage,
    fill = comparison_celltype
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(x * 100),
        "%"
      )
    }
  ) +
  labs(
    title =
      "Breast Cancer Subtype vs Normal Breast Cell Composition",
    x = NULL,
    y = "Cell proportion",
    fill = "Comparison cell type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "celltype_composition_subtype_vs_normal.png"
  ),
  p2,
  width = 11,
  height = 7,
  dpi = 300
)


# ============================================================
# 19. DATASET SUMMARY
# ============================================================

dataset_summary <- bind_rows(
  tumor@meta.data %>%
    count(
      clinical_group = as.character(
        clinical_subtype
      ),
      name = "cells"
    ) %>%
    mutate(
      dataset = "Breast_Cancer"
    ),
  normal@meta.data %>%
    count(
      clinical_group = "Normal",
      name = "cells"
    ) %>%
    mutate(
      dataset = "Normal_Breast"
    )
) %>%
  select(
    dataset,
    clinical_group,
    cells
  )


dataset_summary <- dataset_summary %>%
  group_by(dataset) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  ) %>%
  ungroup()


write.csv(
  dataset_summary,
  file.path(
    output_dir,
    "dataset_and_subtype_cell_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. SAVE COMPARISON-READY OBJECTS
# ============================================================

tumor_output <- file.path(
  output_dir,
  "GSE176078_breast_cancer_comparison_ready_v2.rds"
)

normal_output <- file.path(
  output_dir,
  "GSE113196_normal_breast_comparison_ready_v2.rds"
)


saveRDS(
  tumor,
  tumor_output,
  compress = FALSE
)

saveRDS(
  normal,
  normal_output,
  compress = FALSE
)


# ============================================================
# 21. FINAL VALIDATION
# ============================================================

cat("\n==============================================\n")
cat("FINAL COHORT VALIDATION\n")
cat("==============================================\n\n")

cat(
  "GSE176078 patients: ",
  length(
    unique(
      tumor$patient_id
    )
  ),
  "\n",
  sep = ""
)

cat(
  "GSE176078 cells: ",
  ncol(tumor),
  "\n",
  sep = ""
)

cat(
  "GSE113196 cells: ",
  ncol(normal),
  "\n\n",
  sep = ""
)


cat("Clinical groups:\n")

print(
  table(
    tumor$clinical_subtype
  )
)

cat("\n")


cat("Tumor comparison cell types:\n")

print(
  table(
    tumor$comparison_celltype
  )
)

cat("\n")


cat("Normal comparison cell types:\n")

print(
  table(
    normal$comparison_celltype
  )
)

cat("\n")


# ============================================================
# 22. FINAL OUTPUT SUMMARY
# ============================================================

cat("==============================================\n")
cat("SCRIPT 24 COMPLETE\n")
cat("==============================================\n\n")

cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Key outputs:\n")
cat("1. GSE176078_patient_clinical_subtypes.csv\n")
cat("2. combined_cell_metadata.csv\n")
cat("3. celltype_composition_all_groups.csv\n")
cat("4. celltype_composition_wide.csv\n")
cat("5. celltype_composition_subtype_vs_normal.csv\n")
cat("6. dataset_and_subtype_cell_counts.csv\n")
cat("7. unresolved_tumor_celltypes.csv\n")
cat("8. unresolved_normal_celltypes.csv\n")
cat("9. celltype_composition_all_groups.png\n")
cat("10. celltype_composition_subtype_vs_normal.png\n")
cat("11. GSE176078_breast_cancer_comparison_ready_v2.rds\n")
cat("12. GSE113196_normal_breast_comparison_ready_v2.rds\n\n")

cat("Next step:\n")
cat("Script 25 — subtype-aware disease analysis / DE preparation.\n")
cat("==============================================\n")