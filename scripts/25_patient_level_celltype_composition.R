# ============================================================
# SCRIPT 25 — PATIENT-LEVEL CELL-TYPE COMPOSITION ANALYSIS
# ============================================================
#
# GSE176078:
#   ER+, HER2+, TNBC breast cancer
#
# GSE113196:
#   Normal breast reference
#
# Input:
#   data/processed/GSE176078_breast_cancer_comparison_ready_v2.rds
#   data/processed/GSE113196_normal_breast_comparison_ready_v2.rds
#
# Purpose:
#   1. Calculate patient-level cell-type composition
#   2. Compare ER+, HER2+, and TNBC tumors
#   3. Compare each tumor subtype with normal breast
#   4. Perform patient-level statistical testing
#   5. Generate composition plots
#   6. Produce downstream-ready composition tables
#
# IMPORTANT:
#   Cells are NOT treated as independent biological replicates.
#   Tumor statistics are performed at the patient level.
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat("==============================================\n")
cat("SCRIPT 25 — PATIENT-LEVEL CELL-TYPE COMPOSITION\n")
cat("==============================================\n\n")


# ============================================================
# 1. FILE PATHS
# ============================================================

tumor_file <- "results/cohort_harmonization/GSE176078_breast_cancer_comparison_ready_v2.rds"

normal_file <- "results/cohort_harmonization/GSE113196_normal_breast_comparison_ready_v2.rds"

output_dir <- "results/patient_level_composition"

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
    "Tumor comparison-ready object not found: ",
    tumor_file
  )
}

if (!file.exists(normal_file)) {
  stop(
    "Normal comparison-ready object not found: ",
    normal_file
  )
}


# ============================================================
# 3. LOAD OBJECTS
# ============================================================

cat("Loading tumor object...\n")

tumor <- readRDS(tumor_file)

cat(
  "Tumor cells:",
  ncol(tumor),
  "\n\n"
)


cat("Loading normal breast object...\n")

normal <- readRDS(normal_file)

cat(
  "Normal cells:",
  ncol(normal),
  "\n\n"
)


# ============================================================
# 4. IDENTIFY METADATA COLUMNS
# ============================================================

cat("TUMOR METADATA COLUMNS:\n")
print(colnames(tumor@meta.data))

cat("\nNORMAL METADATA COLUMNS:\n")
print(colnames(normal@meta.data))

cat("\n")


# ============================================================
# 5. FIND REQUIRED METADATA
# ============================================================

# Tumor patient ID

tumor_patient_candidates <- c(
  "patient_id",
  "patient",
  "donor_id",
  "sample_id"
)

tumor_patient_col <- intersect(
  tumor_patient_candidates,
  colnames(tumor@meta.data)
)[1]


# Tumor subtype

tumor_subtype_candidates <- c(
  "clinical_subtype",
  "subtype",
  "clinical_subtype_original"
)

tumor_subtype_col <- intersect(
  tumor_subtype_candidates,
  colnames(tumor@meta.data)
)[1]


# Tumor cell type

tumor_celltype_candidates <- c(
  "celltype_harmonized",
  "celltype_major",
  "celltype_subtype",
  "celltype"
)

tumor_celltype_col <- intersect(
  tumor_celltype_candidates,
  colnames(tumor@meta.data)
)[1]


# Normal individual

normal_patient_candidates <- c(
  "individual",
  "patient_id",
  "donor_id",
  "sample_id"
)

normal_patient_col <- intersect(
  normal_patient_candidates,
  colnames(normal@meta.data)
)[1]


# Normal cell type

normal_celltype_candidates <- c(
  "celltype_harmonized",
  "celltype_v1",
  "celltype",
  "celltype_subtype"
)

normal_celltype_col <- intersect(
  normal_celltype_candidates,
  colnames(normal@meta.data)
)[1]


# ============================================================
# 6. VALIDATE METADATA
# ============================================================

required_columns <- c(
  tumor_patient_col,
  tumor_subtype_col,
  tumor_celltype_col,
  normal_patient_col,
  normal_celltype_col
)

if (any(is.na(required_columns))) {

  stop(
    paste(
      "Could not identify required metadata columns.\n",
      "Tumor patient:",
      tumor_patient_col,
      "\nTumor subtype:",
      tumor_subtype_col,
      "\nTumor cell type:",
      tumor_celltype_col,
      "\nNormal individual:",
      normal_patient_col,
      "\nNormal cell type:",
      normal_celltype_col
    )
  )
}


cat("==============================================\n")
cat("METADATA USED\n")
cat("==============================================\n")

cat(
  "Tumor patient:",
  tumor_patient_col,
  "\n"
)

cat(
  "Tumor subtype:",
  tumor_subtype_col,
  "\n"
)

cat(
  "Tumor cell type:",
  tumor_celltype_col,
  "\n"
)

cat(
  "Normal individual:",
  normal_patient_col,
  "\n"
)

cat(
  "Normal cell type:",
  normal_celltype_col,
  "\n\n"
)


# ============================================================
# 7. BUILD TUMOR CELL METADATA TABLE
# ============================================================

tumor_meta <- tumor@meta.data %>%
  mutate(
    cell_id = rownames(tumor@meta.data),
    patient_id = as.character(.data[[tumor_patient_col]]),
    subtype = as.character(.data[[tumor_subtype_col]]),
    celltype = as.character(.data[[tumor_celltype_col]])
  ) %>%
  select(
    cell_id,
    patient_id,
    subtype,
    celltype
  )


# ============================================================
# 8. BUILD NORMAL CELL METADATA TABLE
# ============================================================

normal_meta <- normal@meta.data %>%
  mutate(
    cell_id = rownames(normal@meta.data),
    patient_id = as.character(.data[[normal_patient_col]]),
    subtype = "Normal",
    celltype = as.character(.data[[normal_celltype_col]])
  ) %>%
  select(
    cell_id,
    patient_id,
    subtype,
    celltype
  )


# ============================================================
# 9. CLEAN LABELS
# ============================================================

tumor_meta <- tumor_meta %>%
  mutate(
    patient_id = trimws(patient_id),
    subtype = trimws(subtype),
    celltype = trimws(celltype)
  )

normal_meta <- normal_meta %>%
  mutate(
    patient_id = trimws(patient_id),
    subtype = "Normal",
    celltype = trimws(celltype)
  )


# ============================================================
# 10. CHECK FOR MISSING VALUES
# ============================================================

cat("Missing tumor patient IDs:",
    sum(is.na(tumor_meta$patient_id) |
          tumor_meta$patient_id == ""),
    "\n"
)

cat("Missing tumor subtypes:",
    sum(is.na(tumor_meta$subtype) |
          tumor_meta$subtype == ""),
    "\n"
)

cat("Missing tumor cell types:",
    sum(is.na(tumor_meta$celltype) |
          tumor_meta$celltype == ""),
    "\n\n"
)


# ============================================================
# 11. PATIENT-LEVEL TUMOR SUBTYPE VALIDATION
# ============================================================

patient_subtype_check <- tumor_meta %>%
  distinct(
    patient_id,
    subtype
  ) %>%
  arrange(patient_id)


cat("==============================================\n")
cat("PATIENT-LEVEL SUBTYPE VALIDATION\n")
cat("==============================================\n")

print(patient_subtype_check)

cat("\n")


patients_with_multiple_subtypes <- patient_subtype_check %>%
  count(patient_id) %>%
  filter(n > 1)


if (nrow(patients_with_multiple_subtypes) > 0) {

  cat(
    "WARNING:",
    nrow(patients_with_multiple_subtypes),
    "patients have multiple subtype labels.\n\n"
  )

  print(patients_with_multiple_subtypes)

  stop(
    "Patient-level clinical subtype inconsistency detected."
  )
}


write.csv(
  patient_subtype_check,
  file.path(
    output_dir,
    "patient_subtype_validation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. COMBINED CELL METADATA
# ============================================================

combined_meta <- bind_rows(
  tumor_meta %>%
    mutate(
      dataset = "GSE176078_Tumor"
    ),
  normal_meta %>%
    mutate(
      dataset = "GSE113196_Normal"
    )
)


write.csv(
  combined_meta,
  file.path(
    output_dir,
    "combined_cell_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. CELL COUNTS PER PATIENT
# ============================================================

patient_cell_counts <- combined_meta %>%
  count(
    dataset,
    subtype,
    patient_id,
    name = "total_cells"
  )


write.csv(
  patient_cell_counts,
  file.path(
    output_dir,
    "patient_total_cell_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. CELL-TYPE COUNTS PER PATIENT
# ============================================================

patient_celltype_counts <- combined_meta %>%
  count(
    dataset,
    subtype,
    patient_id,
    celltype,
    name = "cell_count"
  )


# ============================================================
# 15. CALCULATE PATIENT-LEVEL CELL-TYPE PERCENTAGES
# ============================================================

patient_composition <- patient_celltype_counts %>%
  left_join(
    patient_cell_counts,
    by = c(
      "dataset",
      "subtype",
      "patient_id"
    )
  ) %>%
  mutate(
    percentage = (
      cell_count /
        total_cells
    ) * 100
  )


write.csv(
  patient_celltype_counts,
  file.path(
    output_dir,
    "patient_celltype_counts.csv"
  ),
  row.names = FALSE
)


write.csv(
  patient_composition,
  file.path(
    output_dir,
    "patient_celltype_composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SUMMARY BY CLINICAL GROUP
# ============================================================

group_summary <- patient_composition %>%
  group_by(
    subtype,
    celltype
  ) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    mean_percentage = mean(
      percentage,
      na.rm = TRUE
    ),
    median_percentage = median(
      percentage,
      na.rm = TRUE
    ),
    sd_percentage = sd(
      percentage,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


write.csv(
  group_summary,
  file.path(
    output_dir,
    "group_level_celltype_composition.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. TUMOR-ONLY CELL-TYPE COMPARISON
# ============================================================

tumor_composition <- patient_composition %>%
  filter(
    dataset == "GSE176078_Tumor"
  )


cat("==============================================\n")
cat("TUMOR PATIENT COUNTS\n")
cat("==============================================\n")

print(
  tumor_composition %>%
    distinct(
      subtype,
      patient_id
    ) %>%
    count(subtype)
)

cat("\n")


# ============================================================
# 18. PAIRWISE TUMOR SUBTYPE STATISTICS
# ============================================================

tumor_subtypes <- c(
  "ER_Positive",
  "HER2_Positive",
  "TNBC"
)


pairwise_results <- list()

comparison_number <- 0


for (cell in unique(tumor_composition$celltype)) {

  cell_data <- tumor_composition %>%
    filter(
      celltype == cell
    )

  for (i in seq_len(
    length(tumor_subtypes) - 1
  )) {

    for (j in (i + 1):length(tumor_subtypes)) {

      group1 <- tumor_subtypes[i]
      group2 <- tumor_subtypes[j]

      values1 <- cell_data %>%
        filter(
          subtype == group1
        ) %>%
        pull(percentage)

      values2 <- cell_data %>%
        filter(
          subtype == group2
        ) %>%
        pull(percentage)

      if (
        length(values1) >= 2 &&
        length(values2) >= 2
      ) {

        test <- wilcox.test(
          values1,
          values2,
          exact = FALSE
        )

        comparison_number <-
          comparison_number + 1

        pairwise_results[[comparison_number]] <-
          data.frame(
            celltype = cell,
            comparison = paste(
              group1,
              "vs",
              group2
            ),
            group1 = group1,
            group2 = group2,
            n_group1 = length(values1),
            n_group2 = length(values2),
            median_group1 = median(values1),
            median_group2 = median(values2),
            p_value = test$p.value,
            stringsAsFactors = FALSE
          )
      }
    }
  }
}


tumor_pairwise_results <- bind_rows(
  pairwise_results
)


if (nrow(tumor_pairwise_results) > 0) {

  tumor_pairwise_results <- tumor_pairwise_results %>%
    mutate(
      p_adj = p.adjust(
        p_value,
        method = "BH"
      )
    )
}


write.csv(
  tumor_pairwise_results,
  file.path(
    output_dir,
    "tumor_subtype_pairwise_composition_statistics.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. COMPOSITION PLOT — PATIENT LEVEL
# ============================================================

cat("Generating patient-level composition plot...\n")


plot_data <- patient_composition %>%
  mutate(
    subtype = factor(
      subtype,
      levels = c(
        "Normal",
        "ER_Positive",
        "HER2_Positive",
        "TNBC"
      )
    )
  )


p1 <- ggplot(
  plot_data,
  aes(
    x = subtype,
    y = percentage,
    fill = celltype
  )
) +
  geom_boxplot(
    position = position_dodge(
      width = 0.8
    ),
    outlier.shape = NA
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.5,
    size = 1
  ) +
  labs(
    title = "Patient-level cell-type composition",
    x = "Clinical group",
    y = "Cell-type percentage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )


ggsave(
  file.path(
    output_dir,
    "patient_level_celltype_composition.png"
  ),
  p1,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# 20. COMPOSITION HEATMAP DATA
# ============================================================

composition_wide <- patient_composition %>%
  select(
    patient_id,
    subtype,
    celltype,
    percentage
  ) %>%
  pivot_wider(
    names_from = celltype,
    values_from = percentage,
    values_fill = 0
  )


write.csv(
  composition_wide,
  file.path(
    output_dir,
    "patient_celltype_composition_wide.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. GROUP-LEVEL SUMMARY PLOT
# ============================================================

group_plot_data <- group_summary %>%
  mutate(
    subtype = factor(
      subtype,
      levels = c(
        "Normal",
        "ER_Positive",
        "HER2_Positive",
        "TNBC"
      )
    )
  )


p2 <- ggplot(
  group_plot_data,
  aes(
    x = subtype,
    y = mean_percentage,
    fill = celltype
  )
) +
  geom_col(
    position = "stack"
  ) +
  labs(
    title = "Mean cell-type composition by clinical group",
    x = "Clinical group",
    y = "Mean percentage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )


ggsave(
  file.path(
    output_dir,
    "group_level_celltype_composition.png"
  ),
  p2,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# 22. IDENTIFY MAJOR COMPOSITIONAL DIFFERENCES
# ============================================================

major_differences <- group_summary %>%
  filter(
    !is.na(sd_percentage)
  ) %>%
  arrange(
    subtype,
    desc(mean_percentage)
  )


write.csv(
  major_differences,
  file.path(
    output_dir,
    "celltype_composition_summary_for_interpretation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. SAVE ANALYSIS SUMMARY
# ============================================================

summary_table <- data.frame(
  metric = c(
    "Tumor cells",
    "Normal cells",
    "Tumor patients",
    "Normal individuals",
    "Tumor clinical groups"
  ),
  value = c(
    nrow(tumor_meta),
    nrow(normal_meta),
    n_distinct(
      tumor_meta$patient_id
    ),
    n_distinct(
      normal_meta$patient_id
    ),
    paste(
      sort(
        unique(
          tumor_meta$subtype
        )
      ),
      collapse = ", "
    )
  ),
  stringsAsFactors = FALSE
)


write.csv(
  summary_table,
  file.path(
    output_dir,
    "script25_analysis_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("==============================================\n")
cat("SCRIPT 25 COMPLETE\n")
cat("==============================================\n\n")

cat(
  "Tumor cells:",
  nrow(tumor_meta),
  "\n"
)

cat(
  "Normal cells:",
  nrow(normal_meta),
  "\n"
)

cat(
  "Tumor patients:",
  n_distinct(
    tumor_meta$patient_id
  ),
  "\n"
)

cat(
  "Normal individuals:",
  n_distinct(
    normal_meta$patient_id
  ),
  "\n\n"
)

cat("Clinical groups:\n")

print(
  tumor_meta %>%
    distinct(
      patient_id,
      subtype
    ) %>%
    count(subtype)
)

cat("\n")

cat("Output directory:\n")
cat(
  output_dir,
  "\n\n"
)

cat("Key outputs:\n")

cat(
  "1. patient_subtype_validation.csv\n"
)

cat(
  "2. patient_celltype_counts.csv\n"
)

cat(
  "3. patient_celltype_composition.csv\n"
)

cat(
  "4. group_level_celltype_composition.csv\n"
)

cat(
  "5. tumor_subtype_pairwise_composition_statistics.csv\n"
)

cat(
  "6. patient_celltype_composition_wide.csv\n"
)

cat(
  "7. patient_level_celltype_composition.png\n"
)

cat(
  "8. group_level_celltype_composition.png\n"
)

cat(
  "9. celltype_composition_summary_for_interpretation.csv\n"
)

cat(
  "10. script25_analysis_summary.csv\n"
)

cat("\n==============================================\n")
cat("NEXT STEP: SCRIPT 26 — DISEASE/SUBTYPE EXPRESSION ANALYSIS\n")
cat("==============================================\n")