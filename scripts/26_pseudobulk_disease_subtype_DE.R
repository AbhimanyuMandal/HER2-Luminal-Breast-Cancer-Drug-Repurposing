# ============================================================
# SCRIPT 26 — PSEUDOBULK DISEASE / SUBTYPE DIFFERENTIAL
# EXPRESSION ANALYSIS
#
# GSE176078 breast cancer:
#   ER+
#   HER2+
#   TNBC
#
# GSE113196 normal breast:
#   Normal
#
# Primary analysis:
#   Luminal_Epithelial
#
# Additional comparable compartments:
#   Basal_Epithelial
#   Fibroblast_Stromal
#
# Statistical unit:
#   Patient / individual
#
# Method:
#   Pseudobulk raw counts + edgeR quasi-likelihood model
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(edgeR)
})

cat("==============================================\n")
cat("SCRIPT 26 — PSEUDOBULK DISEASE / SUBTYPE DE\n")
cat("==============================================\n\n")


# ============================================================
# 1. FILE PATHS
# ============================================================

tumor_file <- paste0(
  "results/cohort_harmonization/",
  "GSE176078_breast_cancer_comparison_ready_v2.rds"
)

normal_file <- paste0(
  "results/cohort_harmonization/",
  "GSE113196_normal_breast_comparison_ready_v2.rds"
)

output_dir <- "results/pseudobulk_DE"

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
    "Tumor comparison-ready object not found:\n",
    tumor_file
  )
}

if (!file.exists(normal_file)) {
  stop(
    "Normal comparison-ready object not found:\n",
    normal_file
  )
}


# ============================================================
# 3. CHECK edgeR
# ============================================================

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop(
    "edgeR is required for Script 26.\n",
    "Install it before running this script."
  )
}


# ============================================================
# 4. LOAD OBJECTS
# ============================================================

cat("Loading tumor object...\n")

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


cat("Loading normal breast object...\n")

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
# 5. CHECK RNA COUNTS
# ============================================================

if (!"counts" %in% Layers(tumor[["RNA"]])) {
  stop("Tumor object does not contain an RNA counts layer.")
}

if (!"counts" %in% Layers(normal[["RNA"]])) {
  stop("Normal object does not contain an RNA counts layer.")
}


# ============================================================
# 6. CHECK REQUIRED METADATA
# ============================================================

required_tumor <- c(
  "patient_id",
  "clinical_subtype",
  "comparison_celltype"
)

required_normal <- c(
  "individual",
  "comparison_celltype"
)

missing_tumor <- setdiff(
  required_tumor,
  colnames(tumor@meta.data)
)

missing_normal <- setdiff(
  required_normal,
  colnames(normal@meta.data)
)

if (length(missing_tumor) > 0) {
  stop(
    "Missing tumor metadata:\n",
    paste(missing_tumor, collapse = ", ")
  )
}

if (length(missing_normal) > 0) {
  stop(
    "Missing normal metadata:\n",
    paste(missing_normal, collapse = ", ")
  )
}


# ============================================================
# 7. DEFINE ANALYSIS CELL TYPES
# ============================================================

analysis_celltypes <- c(
  "Luminal_Epithelial",
  "Basal_Epithelial",
  "Fibroblast_Stromal"
)

cat("Analysis cell types:\n")
print(analysis_celltypes)
cat("\n")


# ============================================================
# 8. CHECK CELL-TYPE AVAILABILITY
# ============================================================

cat("==============================================\n")
cat("CELL-TYPE AVAILABILITY\n")
cat("==============================================\n\n")

cat("Tumor:\n")

print(
  table(
    tumor$comparison_celltype,
    useNA = "ifany"
  )
)

cat("\nNormal:\n")

print(
  table(
    normal$comparison_celltype,
    useNA = "ifany"
  )
)

cat("\n")


# ============================================================
# 9. BUILD STANDARDIZED CELL METADATA
# ============================================================

tumor_meta <- tumor@meta.data %>%
  mutate(
    dataset = "Tumor",
    sample_id = as.character(patient_id),
    group = as.character(clinical_subtype),
    celltype = as.character(comparison_celltype)
  ) %>%
  select(
    sample_id,
    group,
    celltype,
    dataset
  )

normal_meta <- normal@meta.data %>%
  mutate(
    dataset = "Normal",
    sample_id = paste0(
      "Normal_",
      as.character(individual)
    ),
    group = "Normal",
    celltype = as.character(comparison_celltype)
  ) %>%
  select(
    sample_id,
    group,
    celltype,
    dataset
  )


# ============================================================
# 10. COMBINE METADATA
# ============================================================

combined_meta <- bind_rows(
  tumor_meta,
  normal_meta
)

cat("Total cells represented:\n")
print(nrow(combined_meta))

cat("\nSamples represented:\n")
print(
  combined_meta %>%
    distinct(sample_id, group, dataset) %>%
    count(dataset, group)
)

cat("\n")


# ============================================================
# 11. PSEUDOBULK FUNCTION
# ============================================================

make_pseudobulk <- function(
  object,
  sample_ids,
  celltype_name,
  sample_column,
  group_column,
  min_cells = 50
) {

  meta <- object@meta.data

  keep <- (
    as.character(meta[[sample_column]]) %in% sample_ids &
    as.character(meta$comparison_celltype) == celltype_name
  )

  cells <- rownames(meta)[keep]

  if (length(cells) == 0) {
    return(NULL)
  }

  sub_meta <- meta[cells, , drop = FALSE]

  counts <- LayerData(
    object,
    assay = "RNA",
    layer = "counts"
  )

  counts <- counts[, cells, drop = FALSE]

  sample_values <- as.character(
    sub_meta[[sample_column]]
  )

  cell_counts <- table(sample_values)

  valid_samples <- names(
    cell_counts[cell_counts >= min_cells]
  )

  if (length(valid_samples) == 0) {
    return(NULL)
  }

  keep_cells <- sample_values %in% valid_samples

  counts <- counts[, keep_cells, drop = FALSE]

  sample_values <- sample_values[keep_cells]

  sub_meta <- sub_meta[
    keep_cells,
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # Efficient sparse pseudobulk aggregation
  #
  # counts:
  #   genes x cells
  #
  # sample_matrix:
  #   cells x samples
  #
  # Matrix multiplication produces:
  #   genes x samples
  #
  # This avoids converting the large sparse count matrix
  # into a dense matrix.
  # ----------------------------------------------------------

  sample_factor <- factor(
    sample_values,
    levels = valid_samples
  )

  sample_matrix <- Matrix::sparse.model.matrix(
    ~ 0 + sample_factor
  )

  colnames(sample_matrix) <- valid_samples

  pseudobulk <- counts %*% sample_matrix

  pseudobulk <- as.matrix(pseudobulk)

  sample_order <- colnames(pseudobulk)

  sample_groups <- sapply(
    sample_order,
    function(id) {

      vals <- unique(
        as.character(
          sub_meta[[group_column]][
            sample_values == id
          ]
        )
      )

      vals[1]
    }
  )

  sample_info <- data.frame(
    sample_id = sample_order,
    group = sample_groups,
    celltype = celltype_name,
    cell_count = as.integer(
      cell_counts[sample_order]
    ),
    stringsAsFactors = FALSE
  )

  rownames(sample_info) <- sample_info$sample_id

  list(
    counts = pseudobulk,
    metadata = sample_info
  )
}


# ============================================================
# 12. EDGE-R DE FUNCTION
# ============================================================

run_edgeR_DE <- function(
  counts,
  sample_metadata,
  group1,
  group2,
  celltype_name
) {

  keep_samples <- sample_metadata$group %in% c(
    group1,
    group2
  )

  sample_metadata <- sample_metadata[
    keep_samples,
    ,
    drop = FALSE
  ]

  counts <- counts[
    ,
    sample_metadata$sample_id,
    drop = FALSE
  ]

  sample_metadata$group <- factor(
    sample_metadata$group,
    levels = c(group2, group1)
  )

  if (length(unique(sample_metadata$group)) < 2) {
    return(NULL)
  }

  if (
    sum(sample_metadata$group == group1) < 2 ||
    sum(sample_metadata$group == group2) < 2
  ) {

    warning(
      "Insufficient biological replicates for ",
      group1,
      " vs ",
      group2,
      " in ",
      celltype_name
    )

    return(NULL)
  }

  dge <- DGEList(
    counts = counts,
    group = sample_metadata$group
  )

  # Filter lowly expressed genes.
  keep_genes <- filterByExpr(
    dge,
    group = sample_metadata$group
  )

  dge <- dge[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  if (nrow(dge) < 100) {

    warning(
      "Too few genes after filtering for ",
      celltype_name,
      ": ",
      nrow(dge)
    )

    return(NULL)
  }

  dge <- calcNormFactors(dge)

  design <- model.matrix(
    ~ group,
    data = sample_metadata
  )

  dge <- estimateDisp(
    dge,
    design
  )

  fit <- glmQLFit(
    dge,
    design,
    robust = TRUE
  )

  qlf <- glmQLFTest(
    fit,
    coef = 2
  )

  result <- topTags(
    qlf,
    n = Inf
  )$table

  result$gene <- rownames(result)

  result$celltype <- celltype_name

  result$group1 <- group1

  result$group2 <- group2

  result$comparison <- paste(
    group1,
    "vs",
    group2
  )

  result$direction <- ifelse(
    result$logFC > 0,
    paste0(group1, "_higher"),
    paste0(group2, "_higher")
  )

  result <- result %>%
    select(
      gene,
      celltype,
      group1,
      group2,
      comparison,
      logFC,
      logCPM,
      F,
      PValue,
      FDR,
      direction,
      everything()
    )

  result
}


# ============================================================
# 13. PREPARE SAMPLE IDS
# ============================================================

tumor_samples <- unique(
  as.character(tumor$patient_id)
)

normal_samples <- unique(
  as.character(normal$individual)
)

cat("Tumor patients:", length(tumor_samples), "\n")
cat("Normal individuals:", length(normal_samples), "\n\n")


# ============================================================
# 14. CREATE PSEUDOBULK DATA
# ============================================================

pseudobulk_objects <- list()

for (ct in analysis_celltypes) {

  cat("==============================================\n")
  cat("PSEUDOBULK:", ct, "\n")
  cat("==============================================\n")

  tumor_pb <- make_pseudobulk(
    object = tumor,
    sample_ids = tumor_samples,
    celltype_name = ct,
    sample_column = "patient_id",
    group_column = "clinical_subtype",
    min_cells = 50
  )

  normal_pb <- make_pseudobulk(
    object = normal,
    sample_ids = normal_samples,
    celltype_name = ct,
    sample_column = "individual",
    group_column = "clinical_subtype",
    min_cells = 50
  )

  if (is.null(tumor_pb)) {
    cat("No usable tumor pseudobulk samples.\n\n")
    next
  }

  if (is.null(normal_pb)) {
    cat("No usable normal pseudobulk samples.\n\n")
    next
  }

  common_genes <- intersect(
    rownames(tumor_pb$counts),
    rownames(normal_pb$counts)
  )

  tumor_counts <- tumor_pb$counts[
    common_genes,
    ,
    drop = FALSE
  ]

  normal_counts <- normal_pb$counts[
    common_genes,
    ,
    drop = FALSE
  ]

  combined_counts <- cbind(
    tumor_counts,
    normal_counts
  )

  combined_metadata <- bind_rows(
    tumor_pb$metadata,
    normal_pb$metadata
  )

  pseudobulk_objects[[ct]] <- list(
    counts = combined_counts,
    metadata = combined_metadata
  )

  cat(
    "Tumor samples:",
    ncol(tumor_counts),
    "\n"
  )

  cat(
    "Normal samples:",
    ncol(normal_counts),
    "\n"
  )

  cat(
    "Common genes:",
    length(common_genes),
    "\n\n"
  )

  saveRDS(
    pseudobulk_objects[[ct]],
    file.path(
      output_dir,
      paste0(
        "pseudobulk_",
        ct,
        ".rds"
      )
    )
  )
}


# ============================================================
# 15. RUN DIFFERENTIAL EXPRESSION
# ============================================================

comparisons <- list(
  c("ER_Positive", "Normal"),
  c("HER2_Positive", "Normal"),
  c("TNBC", "Normal"),
  c("HER2_Positive", "ER_Positive"),
  c("TNBC", "ER_Positive"),
  c("TNBC", "HER2_Positive")
)


all_results <- list()

result_index <- 1


for (ct in names(pseudobulk_objects)) {

  cat("\n")
  cat("##############################################\n")
  cat("CELL TYPE:", ct, "\n")
  cat("##############################################\n")

  pb <- pseudobulk_objects[[ct]]

  for (comparison in comparisons) {

    group1 <- comparison[1]
    group2 <- comparison[2]

    cat(
      "\nRunning:",
      group1,
      "vs",
      group2,
      "\n"
    )

    result <- tryCatch(

      run_edgeR_DE(
        counts = pb$counts,
        sample_metadata = pb$metadata,
        group1 = group1,
        group2 = group2,
        celltype_name = ct
      ),

      error = function(e) {

        warning(
          "DE failed for ",
          ct,
          ": ",
          group1,
          " vs ",
          group2,
          "\n",
          e$message
        )

        NULL
      }
    )

    if (!is.null(result)) {

      all_results[[result_index]] <- result

      result_index <- result_index + 1

      safe_name <- gsub(
        "[^A-Za-z0-9]+",
        "_",
        paste(
          ct,
          group1,
          "vs",
          group2
        )
      )

      write.csv(
        result,
        file.path(
          output_dir,
          paste0(
            safe_name,
            ".csv"
          )
        ),
        row.names = FALSE
      )

      # Significant genes

      significant <- result %>%
        filter(
          FDR < 0.05,
          abs(logFC) >= 1
        )

      write.csv(
        significant,
        file.path(
          output_dir,
          paste0(
            safe_name,
            "_significant.csv"
          )
        ),
        row.names = FALSE
      )

      cat(
        "Significant genes:",
        nrow(significant),
        "\n"
      )
    }
  }
}


# ============================================================
# 16. COMBINE ALL DE RESULTS
# ============================================================

if (length(all_results) > 0) {

  combined_DE <- bind_rows(
    all_results
  )

  write.csv(
    combined_DE,
    file.path(
      output_dir,
      "all_pseudobulk_DE_results.csv"
    ),
    row.names = FALSE
  )

  significant_all <- combined_DE %>%
    filter(
      FDR < 0.05,
      abs(logFC) >= 1
    )

  write.csv(
    significant_all,
    file.path(
      output_dir,
      "all_significant_DE_genes.csv"
    ),
    row.names = FALSE
  )

} else {

  stop(
    "No differential-expression results were generated."
  )
}


# ============================================================
# 17. DE SUMMARY
# ============================================================

de_summary <- combined_DE %>%
  group_by(
    celltype,
    comparison
  ) %>%
  summarise(
    genes_tested = n(),
    significant_FDR05 = sum(
      FDR < 0.05,
      na.rm = TRUE
    ),
    significant_FDR05_logFC1 = sum(
      FDR < 0.05 &
        abs(logFC) >= 1,
      na.rm = TRUE
    ),
    upregulated = sum(
      FDR < 0.05 &
        logFC >= 1,
      na.rm = TRUE
    ),
    downregulated = sum(
      FDR < 0.05 &
        logFC <= -1,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  de_summary,
  file.path(
    output_dir,
    "DE_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. TOP GENES
# ============================================================

top_genes <- combined_DE %>%
  filter(
    FDR < 0.05
  ) %>%
  group_by(
    celltype,
    comparison
  ) %>%
  arrange(
    FDR,
    desc(abs(logFC))
  ) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top_genes,
  file.path(
    output_dir,
    "top20_DE_genes_per_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. SIMPLE VOLCANO PLOTS
# ============================================================

plot_volcano <- function(
  result,
  title,
  filename
) {

  result <- result %>%
    mutate(
      significance = case_when(
        FDR < 0.05 &
          logFC >= 1 ~ "Up",
        FDR < 0.05 &
          logFC <= -1 ~ "Down",
        TRUE ~ "Not significant"
      ),
      neg_log10_FDR = -log10(
        pmax(FDR, 1e-300)
      )
    )

  p <- ggplot(
    result,
    aes(
      x = logFC,
      y = neg_log10_FDR
    )
  ) +
    geom_point(
      aes(
        alpha = significance
      ),
      size = 1
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(FDR)",
      alpha = NULL
    ) +
    theme_minimal()

  ggsave(
    filename,
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# Generate plots only for the primary luminal compartment.

if ("Luminal_Epithelial" %in% names(pseudobulk_objects)) {

  luminal_results <- combined_DE %>%
    filter(
      celltype == "Luminal_Epithelial"
    )

  unique_comparisons <- unique(
    luminal_results$comparison
  )

  for (cmp in unique_comparisons) {

    tmp <- luminal_results %>%
      filter(
        comparison == cmp
      )

    safe_name <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      cmp
    )

    plot_volcano(
      tmp,
      paste(
        "Luminal epithelial:",
        cmp
      ),
      file.path(
        output_dir,
        paste0(
          "volcano_luminal_",
          safe_name,
          ".png"
        )
      )
    )
  }
}


# ============================================================
# 20. FINAL SUMMARY
# ============================================================

cat("\n")
cat("==============================================\n")
cat("SCRIPT 26 COMPLETE\n")
cat("==============================================\n\n")

cat(
  "Tumor patients:",
  length(tumor_samples),
  "\n"
)

cat(
  "Normal individuals:",
  length(normal_samples),
  "\n\n"
)

cat("Cell types analyzed:\n")
print(names(pseudobulk_objects))

cat("\nDE summary:\n")
print(de_summary)

cat("\nOutput directory:\n")
cat(
  output_dir,
  "\n\n"
)

cat("Important interpretation rule:\n")
cat(
  "Patient/individual is the biological replicate.\n"
)

cat(
  "HER2-positive comparisons contain only 5 tumor patients,\n"
  ,"so these results should be treated cautiously.\n"
)

cat(
  "\nPrimary drug-repurposing compartment:\n"
)

cat(
  "Luminal_Epithelial\n"
)

cat(
  "\nNext step after reviewing Script 26:\n"
)

cat(
  "Pathway/signature analysis followed by disease-signature\n"
  ,"construction for drug repurposing.\n"
)