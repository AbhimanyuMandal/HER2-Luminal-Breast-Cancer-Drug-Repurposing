# ============================================================
# SCRIPT 27 — HER2+ LUMINAL DISEASE SIGNATURE QC
# ============================================================
#
# Purpose:
#   1. Load HER2+ Luminal_Epithelial vs Normal DE results
#   2. Preserve the complete DE result
#   3. Apply statistical/effect-size thresholds
#   4. Identify immune-associated genes
#   5. Identify epithelial/luminal markers
#   6. Construct broad and strict disease signatures
#   7. Rank genes for downstream drug repurposing
#   8. Generate QC summaries and plots
#
# Primary comparison:
#   HER2_Positive Luminal_Epithelial
#       vs
#   Normal Luminal_Epithelial
#
# IMPORTANT:
#   This script does NOT redo differential expression.
#   This script does NOT alter the original DE results.
#
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

cat("==============================================\n")
cat("SCRIPT 27 — HER2+ LUMINAL DISEASE SIGNATURE QC\n")
cat("==============================================\n\n")


# ============================================================
# 1. INPUT / OUTPUT PATHS
# ============================================================

input_file <- paste0(
  "results/pseudobulk_DE/",
  "Luminal_Epithelial_HER2_Positive_vs_Normal.csv"
)

output_dir <- "results/disease_signatures/HER2_Luminal"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "HER2+ vs Normal DE file not found:\n",
    input_file
  )
}


# ============================================================
# 2. LOAD DE RESULTS
# ============================================================

cat("Loading HER2+ Luminal DE results...\n")

de <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

cat(
  "Genes loaded:",
  nrow(de),
  "\n\n"
)


# ============================================================
# 3. CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "gene",
  "logFC",
  "FDR",
  "PValue"
)

missing_columns <- setdiff(
  required_columns,
  colnames(de)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ============================================================
# 4. BASIC DE QC
# ============================================================

cat("==============================================\n")
cat("BASIC DE QC\n")
cat("==============================================\n\n")

cat(
  "Total genes tested:",
  nrow(de),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(
    de$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FDR < 0.05 & |logFC| >= 1:",
  sum(
    de$FDR < 0.05 &
      abs(de$logFC) >= 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Upregulated:",
  sum(
    de$FDR < 0.05 &
      de$logFC >= 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Downregulated:",
  sum(
    de$FDR < 0.05 &
      de$logFC <= -1,
    na.rm = TRUE
  ),
  "\n\n"
)


# ============================================================
# 5. DEFINE SIGNIFICANT DE GENES
# ============================================================

sig <- de %>%
  filter(
    !is.na(FDR),
    !is.na(logFC),
    FDR < 0.05,
    abs(logFC) >= 1
  )

cat(
  "Significant genes for signature construction:",
  nrow(sig),
  "\n\n"
)


# ============================================================
# 6. IMMUNE-ASSOCIATED GENE LIST
# ============================================================
#
# These genes are NOT deleted from the original DE results.
#
# They are flagged because our earlier QC showed a small
# immune signal within the HER2+ luminal compartment.
#
# ============================================================

immune_markers <- c(
  "PTPRC",
  "CD3D",
  "CD3E",
  "CD3G",
  "CD7",
  "CD8A",
  "CD8B",
  "IL7R",
  "LST1",
  "TYROBP",
  "FCER1G",
  "LYZ",
  "TRBC1",
  "TRBC2",
  "NKG7",
  "GNLY",
  "MS4A1",
  "CD79A",
  "CD37",
  "LAPTM5",
  "CD14",
  "CTSS",
  "HLA-DRA",
  "HLA-DPB1",
  "HLA-DPA1",
  "HLA-DQB1",
  "CST3",
  "S100A8",
  "S100A9"
)

immune_found <- intersect(
  immune_markers,
  sig$gene
)

cat("==============================================\n")
cat("IMMUNE SIGNAL QC\n")
cat("==============================================\n\n")

cat(
  "Significant immune-associated genes:",
  length(immune_found),
  "\n\n"
)

if (length(immune_found) > 0) {

  immune_table <- sig %>%
    filter(
      gene %in% immune_found
    ) %>%
    arrange(
      desc(abs(logFC))
    )

  print(
    immune_table[
      ,
      c(
        "gene",
        "logFC",
        "FDR",
        "PValue"
      )
    ],
    row.names = FALSE
  )

} else {

  cat(
    "No immune-associated genes detected.\n"
  )
}

cat("\n")


# ============================================================
# 7. EPITHELIAL / LUMINAL QC
# ============================================================

luminal_markers <- c(
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "KRT7",
  "MUC1",
  "KRT20",
  "KRT23",
  "KRT24",
  "KRT15"
)

basal_markers <- c(
  "KRT5",
  "KRT6A",
  "KRT6B",
  "KRT14",
  "KRT17"
)

epithelial_markers <- unique(
  c(
    luminal_markers,
    basal_markers
  )
)

epithelial_found <- intersect(
  epithelial_markers,
  sig$gene
)

cat("==============================================\n")
cat("EPITHELIAL / LUMINAL QC\n")
cat("==============================================\n\n")

cat(
  "Significant epithelial markers:",
  length(epithelial_found),
  "\n\n"
)

if (length(epithelial_found) > 0) {

  epithelial_table <- sig %>%
    filter(
      gene %in% epithelial_found
    ) %>%
    arrange(
      desc(abs(logFC))
    )

  print(
    epithelial_table[
      ,
      c(
        "gene",
        "logFC",
        "FDR",
        "PValue"
      )
    ],
    row.names = FALSE
  )

} else {

  cat(
    "No canonical epithelial markers reached the threshold.\n"
  )
}

cat("\n")


# ============================================================
# 8. FLAG GENES
# ============================================================

sig <- sig %>%
  mutate(
    immune_flag = gene %in% immune_markers,
    luminal_flag = gene %in% luminal_markers,
    epithelial_flag = gene %in% epithelial_markers
  )


# ============================================================
# 9. BROAD SIGNATURE
# ============================================================
#
# Broad signature:
#   FDR < 0.05
#   |logFC| >= 1
#
# This preserves all statistically significant biology.
#
# ============================================================

broad_signature <- sig %>%
  arrange(
    desc(abs(logFC))
  )

broad_up <- broad_signature %>%
  filter(
    logFC >= 1
  ) %>%
  arrange(
    desc(logFC)
  )

broad_down <- broad_signature %>%
  filter(
    logFC <= -1
  ) %>%
  arrange(
    logFC
  )


# ============================================================
# 10. STRICT DRUG-REPURPOSING SIGNATURE
# ============================================================
#
# Strict signature excludes genes that are obvious immune
# markers.
#
# We retain epithelial/basal/luminal genes because these may
# represent real tumor-state biology.
#
# ============================================================

strict_signature <- sig %>%
  filter(
    !immune_flag
  ) %>%
  arrange(
    desc(abs(logFC))
  )

strict_up <- strict_signature %>%
  filter(
    logFC >= 1
  ) %>%
  arrange(
    desc(logFC)
  )

strict_down <- strict_signature %>%
  filter(
    logFC <= -1
  ) %>%
  arrange(
    logFC
  )


# ============================================================
# 11. ADD RANKING
# ============================================================

broad_up <- broad_up %>%
  mutate(
    rank = row_number(),
    direction = "UP"
  )

broad_down <- broad_down %>%
  mutate(
    rank = row_number(),
    direction = "DOWN"
  )

strict_up <- strict_up %>%
  mutate(
    rank = row_number(),
    direction = "UP"
  )

strict_down <- strict_down %>%
  mutate(
    rank = row_number(),
    direction = "DOWN"
  )


# ============================================================
# 12. SAVE ORIGINAL SIGNIFICANT DE
# ============================================================

write.csv(
  sig,
  file.path(
    output_dir,
    "HER2_Luminal_significant_DE_all.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. SAVE BROAD SIGNATURE
# ============================================================

write.csv(
  broad_signature,
  file.path(
    output_dir,
    "HER2_Luminal_broad_signature.csv"
  ),
  row.names = FALSE
)

write.csv(
  broad_up,
  file.path(
    output_dir,
    "HER2_Luminal_broad_UP.csv"
  ),
  row.names = FALSE
)

write.csv(
  broad_down,
  file.path(
    output_dir,
    "HER2_Luminal_broad_DOWN.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SAVE STRICT SIGNATURE
# ============================================================

write.csv(
  strict_signature,
  file.path(
    output_dir,
    "HER2_Luminal_strict_signature.csv"
  ),
  row.names = FALSE
)

write.csv(
  strict_up,
  file.path(
    output_dir,
    "HER2_Luminal_strict_UP.csv"
  ),
  row.names = FALSE
)

write.csv(
  strict_down,
  file.path(
    output_dir,
    "HER2_Luminal_strict_DOWN.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. IMMUNE EXCLUSION TABLE
# ============================================================

immune_excluded <- sig %>%
  filter(
    immune_flag
  ) %>%
  arrange(
    desc(abs(logFC))
  )

write.csv(
  immune_excluded,
  file.path(
    output_dir,
    "HER2_Luminal_immune_genes_excluded_from_strict_signature.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SIGNATURE SUMMARY
# ============================================================

summary_table <- data.frame(
  metric = c(
    "Genes_tested",
    "FDR_lt_0.05",
    "FDR_lt_0.05_abs_logFC_ge_1",
    "Broad_UP",
    "Broad_DOWN",
    "Broad_total",
    "Immune_genes_flagged",
    "Strict_UP",
    "Strict_DOWN",
    "Strict_total"
  ),
  value = c(
    nrow(de),

    sum(
      de$FDR < 0.05,
      na.rm = TRUE
    ),

    nrow(sig),

    nrow(broad_up),

    nrow(broad_down),

    nrow(broad_signature),

    length(immune_found),

    nrow(strict_up),

    nrow(strict_down),

    nrow(strict_signature)
  )
)

write.csv(
  summary_table,
  file.path(
    output_dir,
    "HER2_Luminal_signature_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. TOP STRICT UPREGULATED GENES
# ============================================================

top_up <- head(
  strict_up,
  30
)

write.csv(
  top_up,
  file.path(
    output_dir,
    "HER2_Luminal_top30_UP.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. TOP STRICT DOWNREGULATED GENES
# ============================================================

top_down <- head(
  strict_down,
  30
)

write.csv(
  top_down,
  file.path(
    output_dir,
    "HER2_Luminal_top30_DOWN.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. VOLCANO-STYLE QC PLOT
# ============================================================

plot_data <- de %>%
  mutate(
    significance = case_when(
      FDR < 0.05 &
        logFC >= 1 ~ "Up",
      FDR < 0.05 &
        logFC <= -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

p <- ggplot(
  plot_data,
  aes(
    x = logFC,
    y = -log10(FDR),
    alpha = significance
  )
) +
  geom_point(
    size = 1.2
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
    title = "HER2+ Luminal vs Normal",
    subtitle = "Pseudobulk differential expression",
    x = "logFC",
    y = "-log10(FDR)",
    alpha = "Category"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "HER2_Luminal_DE_volcano_QC.png"
  ),
  p,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 20. FINAL REPORT
# ============================================================

cat("\n")
cat("==============================================\n")
cat("SCRIPT 27 COMPLETE\n")
cat("==============================================\n\n")

cat(
  "Genes tested:",
  nrow(de),
  "\n"
)

cat(
  "Significant genes (FDR < 0.05):",
  sum(
    de$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Broad signature:",
  nrow(broad_signature),
  "genes\n"
)

cat(
  "Broad UP:",
  nrow(broad_up),
  "\n"
)

cat(
  "Broad DOWN:",
  nrow(broad_down),
  "\n"
)

cat(
  "Immune genes flagged:",
  length(immune_found),
  "\n"
)

cat(
  "Strict signature:",
  nrow(strict_signature),
  "genes\n"
)

cat(
  "Strict UP:",
  nrow(strict_up),
  "\n"
)

cat(
  "Strict DOWN:",
  nrow(strict_down),
  "\n\n"
)

cat("Output directory:\n")
cat(
  output_dir,
  "\n\n"
)

cat("Key outputs:\n")
cat("1. HER2_Luminal_significant_DE_all.csv\n")
cat("2. HER2_Luminal_broad_UP.csv\n")
cat("3. HER2_Luminal_broad_DOWN.csv\n")
cat("4. HER2_Luminal_strict_UP.csv\n")
cat("5. HER2_Luminal_strict_DOWN.csv\n")
cat("6. HER2_Luminal_immune_genes_excluded_from_strict_signature.csv\n")
cat("7. HER2_Luminal_signature_summary.csv\n")
cat("8. HER2_Luminal_top30_UP.csv\n")
cat("9. HER2_Luminal_top30_DOWN.csv\n")
cat("10. HER2_Luminal_DE_volcano_QC.png\n\n")

cat("IMPORTANT:\n")
cat(
  "The broad signature preserves all significant DE genes.\n"
)

cat(
  "The strict signature excludes only predefined immune-associated genes.\n"
)

cat(
  "Do NOT use the strict signature for biological interpretation without\n"
)

cat(
  "reviewing the excluded immune genes first.\n"
)

cat("\nNext step: Script 28 — pathway/signature analysis.\n")