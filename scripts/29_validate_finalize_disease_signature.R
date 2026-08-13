# ============================================================
# SCRIPT 29 — VALIDATE + FINALIZE HER2+ LUMINAL DISEASE SIGNATURE
# ============================================================
#
# Purpose:
#
#   1. Validate robustness of the Script 27 signature
#   2. Compare broad vs strict signatures
#   3. Define Core and Extended disease signatures
#   4. Preserve UP / DOWN direction
#   5. Flag immune-associated genes
#   6. Add pathway/program context from Script 28
#   7. Perform contextual HER2/luminal marker QC
#   8. Generate final drug-repurposing-ready signatures
#
# IMPORTANT:
#
#   This script does NOT perform differential expression.
#   This script does NOT select drugs.
#   This script does NOT remove genes based on biological preference.
#
# The final signature will be used by Script 30.
#
# ============================================================


suppressPackageStartupMessages({

  library(dplyr)
  library(ggplot2)

})


cat("==============================================\n")
cat("SCRIPT 29 — DISEASE SIGNATURE VALIDATION\n")
cat("==============================================\n\n")


# ============================================================
# 1. DIRECTORIES
# ============================================================

signature_dir <-
  "results/disease_signatures/HER2_Luminal"

pathway_dir <-
  "results/pathway_analysis/HER2_Luminal"

output_dir <-
  "results/final_disease_signature/HER2_Luminal"


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. INPUT FILES
# ============================================================

strict_up_file <- file.path(
  signature_dir,
  "HER2_Luminal_strict_UP.csv"
)

strict_down_file <- file.path(
  signature_dir,
  "HER2_Luminal_strict_DOWN.csv"
)

broad_up_file <- file.path(
  signature_dir,
  "HER2_Luminal_broad_UP.csv"
)

broad_down_file <- file.path(
  signature_dir,
  "HER2_Luminal_broad_DOWN.csv"
)

all_de_file <- file.path(
  signature_dir,
  "HER2_Luminal_significant_DE_all.csv"
)

immune_file <- file.path(
  signature_dir,
  "HER2_Luminal_immune_genes_excluded_from_strict_signature.csv"
)

go_up_file <- file.path(
  pathway_dir,
  "HER2_Luminal_GO_BP_UP.csv"
)

go_down_file <- file.path(
  pathway_dir,
  "HER2_Luminal_GO_BP_DOWN.csv"
)


required_files <- c(
  strict_up_file,
  strict_down_file,
  broad_up_file,
  broad_down_file,
  all_de_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {

  stop(
    paste0(
      "\nMissing required Script 27 file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      "\n\nRun Script 27 first.\n"
    )
  )

}


# ============================================================
# 3. LOAD DATA
# ============================================================

cat("Loading Script 27 disease-signature results...\n\n")


strict_up <- read.csv(
  strict_up_file,
  stringsAsFactors = FALSE
)

strict_down <- read.csv(
  strict_down_file,
  stringsAsFactors = FALSE
)

broad_up <- read.csv(
  broad_up_file,
  stringsAsFactors = FALSE
)

broad_down <- read.csv(
  broad_down_file,
  stringsAsFactors = FALSE
)

all_de <- read.csv(
  all_de_file,
  stringsAsFactors = FALSE
)


cat(
  "Strict UP:",
  nrow(strict_up),
  "\n"
)

cat(
  "Strict DOWN:",
  nrow(strict_down),
  "\n"
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
  "Significant DE:",
  nrow(all_de),
  "\n\n"
)


# ============================================================
# 4. CLEAN GENE LISTS
# ============================================================

strict_up_genes <- unique(
  strict_up$gene[
    !is.na(strict_up$gene) &
      strict_up$gene != ""
  ]
)

strict_down_genes <- unique(
  strict_down$gene[
    !is.na(strict_down$gene) &
      strict_down$gene != ""
  ]
)

broad_up_genes <- unique(
  broad_up$gene[
    !is.na(broad_up$gene) &
      broad_up$gene != ""
  ]
)

broad_down_genes <- unique(
  broad_down$gene[
    !is.na(broad_down$gene) &
      broad_down$gene != ""
  ]
)


# ============================================================
# 5. ROBUSTNESS CHECK
# ============================================================

cat("==============================================\n")
cat("BROAD vs STRICT ROBUSTNESS\n")
cat("==============================================\n\n")


up_overlap <- intersect(
  broad_up_genes,
  strict_up_genes
)

down_overlap <- intersect(
  broad_down_genes,
  strict_down_genes
)


up_retention <-
  100 *
  length(up_overlap) /
  length(broad_up_genes)


down_retention <-
  100 *
  length(down_overlap) /
  length(broad_down_genes)


robustness <- data.frame(

  direction = c(
    "UP",
    "DOWN"
  ),

  broad_genes = c(
    length(broad_up_genes),
    length(broad_down_genes)
  ),

  strict_genes = c(
    length(strict_up_genes),
    length(strict_down_genes)
  ),

  overlapping_genes = c(
    length(up_overlap),
    length(down_overlap)
  ),

  retention_percent = c(
    up_retention,
    down_retention
  )

)


print(
  robustness,
  row.names = FALSE
)


write.csv(
  robustness,
  file.path(
    output_dir,
    "signature_robustness.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. CREATE EXTENDED SIGNATURE
# ============================================================
#
# Extended = all strict significant genes.
#
# This is the broadest defensible disease signature after
# the predefined immune-gene exclusion performed in Script 27.
#
# ============================================================

extended_up <- strict_up %>%
  mutate(
    direction = "UP"
  )

extended_down <- strict_down %>%
  mutate(
    direction = "DOWN"
  )


extended_signature <- bind_rows(
  extended_up,
  extended_down
)


# ============================================================
# 7. CREATE CORE SIGNATURE
# ============================================================
#
# Core signature:
#
#   FDR <= 0.01
#   AND
#   |logFC| >= 1
#
# This is intentionally stricter than the Script 27 threshold.
#
# It is NOT a replacement for the Extended signature.
#
# Both will be carried forward to Script 30.
#
# ============================================================

core_signature <- extended_signature %>%
  filter(
    FDR <= 0.01,
    abs(logFC) >= 1
  )


core_up <- core_signature %>%
  filter(
    direction == "UP"
  )

core_down <- core_signature %>%
  filter(
    direction == "DOWN"
  )


cat("\n==============================================\n")
cat("FINAL SIGNATURE SIZES\n")
cat("==============================================\n\n")


cat(
  "Extended UP:",
  nrow(extended_up),
  "\n"
)

cat(
  "Extended DOWN:",
  nrow(extended_down),
  "\n"
)

cat(
  "Extended TOTAL:",
  nrow(extended_signature),
  "\n\n"
)

cat(
  "Core UP:",
  nrow(core_up),
  "\n"
)

cat(
  "Core DOWN:",
  nrow(core_down),
  "\n"
)

cat(
  "Core TOTAL:",
  nrow(core_signature),
  "\n\n"
)


# ============================================================
# 8. GENE-LEVEL RANKING SCORE
# ============================================================
#
# Ranking score:
#
#   sign(logFC) × -log10(PValue)
#
# Larger positive values:
#   strongly UP
#
# More negative values:
#   strongly DOWN
#
# This ranking is useful for connectivity-based drug
# repurposing later.
#
# ============================================================

ranked_signature <- extended_signature %>%

  mutate(

    rank_score =
      sign(logFC) *
      (-log10(
        pmax(
          PValue,
          .Machine$double.xmin
        )
      ))

  ) %>%

  arrange(
    desc(rank_score)
  )


# ============================================================
# 9. FLAG IMMUNE-ASSOCIATED GENES
# ============================================================
#
# We DO NOT remove these again.
#
# Script 27 already generated the strict signature.
#
# Here we simply document them.
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
  "LAPTM5"

)


ranked_signature$immune_flag <-
  ranked_signature$gene %in%
  immune_markers


cat(
  "Immune-associated genes retained in extended signature:",
  sum(
    ranked_signature$immune_flag
  ),
  "\n\n"
)


# ============================================================
# 10. CONTEXTUAL HER2 / LUMINAL MARKER CHECK
# ============================================================
#
# These genes are NOT used to filter the signature.
#
# They are simply reported for interpretation.
#
# ============================================================

context_markers <- c(

  # HER2 / ERBB2 biology
  "ERBB2",
  "GRB7",
  "EGFR",

  # Luminal / epithelial identity
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "KRT7",
  "MUC1",

  # Hormone-related
  "ESR1",
  "PGR",
  "AR",

  # Luminal transcriptional regulators
  "GATA3",
  "FOXA1",

  # Proliferation
  "MKI67",
  "TOP2A",
  "UBE2C"

)


context_results <- all_de %>%

  filter(
    gene %in% context_markers
  ) %>%

  select(
    gene,
    logFC,
    FDR,
    PValue
  ) %>%

  mutate(
    marker_category = case_when(

      gene %in% c(
        "ERBB2",
        "GRB7",
        "EGFR"
      ) ~ "HER2_signaling",

      gene %in% c(
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT7",
        "MUC1"
      ) ~ "Epithelial_Luminal",

      gene %in% c(
        "ESR1",
        "PGR",
        "AR"
      ) ~ "Hormone_Receptor",

      gene %in% c(
        "GATA3",
        "FOXA1"
      ) ~ "Luminal_Transcriptional",

      gene %in% c(
        "MKI67",
        "TOP2A",
        "UBE2C"
      ) ~ "Proliferation",

      TRUE ~ "Other"

    )
  )


write.csv(
  context_results,
  file.path(
    output_dir,
    "contextual_HER2_luminal_marker_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. PATHWAY PROGRAM ANNOTATION
# ============================================================
#
# If Script 28 GO files are available, extract the genes
# contributing to significant GO pathways.
#
# This provides pathway context for the final signature.
#
# ============================================================

parse_geneID <- function(
  gene_string
) {

  if (
    is.na(gene_string) ||
    gene_string == ""
  ) {

    return(
      character()
    )

  }

  unique(
    unlist(
      strsplit(
        gene_string,
        "/"
      )
    )
  )

}


pathway_gene_table <- data.frame()


if (
  file.exists(go_up_file) &&
  file.exists(go_down_file)
) {

  cat(
    "Loading significant GO pathway gene memberships...\n"
  )


  go_up <- read.csv(
    go_up_file,
    stringsAsFactors = FALSE
  )

  go_down <- read.csv(
    go_down_file,
    stringsAsFactors = FALSE
  )


  if (
    nrow(go_up) > 0 &&
    "geneID" %in% colnames(go_up)
  ) {

    go_up_rows <- lapply(
      seq_len(nrow(go_up)),
      function(i) {

        genes <- parse_geneID(
          go_up$geneID[i]
        )

        if (length(genes) == 0) {
          return(NULL)
        }

        data.frame(

          gene = genes,

          pathway =
            go_up$Description[i],

          direction =
            "UP",

          stringsAsFactors = FALSE

        )

      }
    )

    go_up_table <- bind_rows(
      go_up_rows
    )

  } else {

    go_up_table <- data.frame()

  }


  if (
    nrow(go_down) > 0 &&
    "geneID" %in% colnames(go_down)
  ) {

    go_down_rows <- lapply(
      seq_len(nrow(go_down)),
      function(i) {

        genes <- parse_geneID(
          go_down$geneID[i]
        )

        if (length(genes) == 0) {
          return(NULL)
        }

        data.frame(

          gene = genes,

          pathway =
            go_down$Description[i],

          direction =
            "DOWN",

          stringsAsFactors = FALSE

        )

      }
    )

    go_down_table <- bind_rows(
      go_down_rows
    )

  } else {

    go_down_table <- data.frame()

  }


  pathway_gene_table <- bind_rows(
    go_up_table,
    go_down_table
  )


  if (
    nrow(pathway_gene_table) > 0
  ) {

    pathway_gene_table <-
      pathway_gene_table %>%

      filter(
        gene %in%
          extended_signature$gene
      ) %>%

      distinct()

  }

}


# ============================================================
# 12. GENE-LEVEL PATHWAY SUMMARY
# ============================================================

if (
  nrow(pathway_gene_table) > 0
) {

  pathway_summary <- pathway_gene_table %>%

    group_by(
      gene
    ) %>%

    summarise(

      GO_pathways =
        n_distinct(
          pathway
        ),

      GO_UP_pathways =
        n_distinct(
          pathway[
            direction == "UP"
          ]
        ),

      GO_DOWN_pathways =
        n_distinct(
          pathway[
            direction == "DOWN"
          ]
        ),

      pathway_names =
        paste(
          unique(pathway),
          collapse = " | "
        ),

      .groups = "drop"

    )


  ranked_signature <-
    ranked_signature %>%

    left_join(
      pathway_summary,
      by = "gene"
    )

} else {

  ranked_signature$GO_pathways <- 0
  ranked_signature$GO_UP_pathways <- 0
  ranked_signature$GO_DOWN_pathways <- 0
  ranked_signature$pathway_names <- NA_character_

}


ranked_signature$GO_pathways[
  is.na(
    ranked_signature$GO_pathways
  )
] <- 0

ranked_signature$GO_UP_pathways[
  is.na(
    ranked_signature$GO_UP_pathways
  )
] <- 0

ranked_signature$GO_DOWN_pathways[
  is.na(
    ranked_signature$GO_DOWN_pathways
  )
] <- 0


# ============================================================
# 13. SAVE FINAL SIGNATURES
# ============================================================

write.csv(
  extended_up,
  file.path(
    output_dir,
    "HER2_Luminal_EXTENDED_UP.csv"
  ),
  row.names = FALSE
)

write.csv(
  extended_down,
  file.path(
    output_dir,
    "HER2_Luminal_EXTENDED_DOWN.csv"
  ),
  row.names = FALSE
)

write.csv(
  extended_signature,
  file.path(
    output_dir,
    "HER2_Luminal_EXTENDED_SIGNATURE.csv"
  ),
  row.names = FALSE
)


write.csv(
  core_up,
  file.path(
    output_dir,
    "HER2_Luminal_CORE_UP.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_down,
  file.path(
    output_dir,
    "HER2_Luminal_CORE_DOWN.csv"
  ),
  row.names = FALSE
)

write.csv(
  core_signature,
  file.path(
    output_dir,
    "HER2_Luminal_CORE_SIGNATURE.csv"
  ),
  row.names = FALSE
)


write.csv(
  ranked_signature,
  file.path(
    output_dir,
    "HER2_Luminal_FINAL_RANKED_SIGNATURE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SAVE PATHWAY-GENE TABLE
# ============================================================

if (
  nrow(pathway_gene_table) > 0
) {

  write.csv(
    pathway_gene_table,
    file.path(
      output_dir,
      "HER2_Luminal_signature_GO_pathway_membership.csv"
    ),
    row.names = FALSE
  )

}


# ============================================================
# 15. SIGNATURE QC PLOT
# ============================================================

plot_data <- all_de %>%

  mutate(

    significance =
      case_when(

        FDR <= 0.01 &
          abs(logFC) >= 1
        ~ "Core",

        FDR < 0.05 &
          abs(logFC) >= 1
        ~ "Extended",

        TRUE
        ~ "Other"

      )

  )


p <- ggplot(
  plot_data,
  aes(
    x = logFC,
    y = -log10(
      pmax(
        FDR,
        .Machine$double.xmin
      )
    ),
    alpha = significance
  )
) +

  geom_point(
    size = 1.4
  ) +

  geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept =
      -log10(0.05),
    linetype = "dashed"
  ) +

  scale_alpha_manual(
    values = c(
      "Core" = 1,
      "Extended" = 0.6,
      "Other" = 0.2
    )
  ) +

  labs(
    title =
      "HER2+ Luminal Disease Signature QC",

    subtitle =
      "Core: FDR ≤ 0.01 and |logFC| ≥ 1",

    x = "logFC",

    y = "-log10(FDR)",

    alpha = "Signature class"
  ) +

  theme_bw()


ggsave(
  file.path(
    output_dir,
    "HER2_Luminal_final_signature_QC.png"
  ),
  p,
  width = 10,
  height = 8,
  dpi = 300
)


# ============================================================
# 16. FINAL SUMMARY TABLE
# ============================================================

final_summary <- data.frame(

  metric = c(

    "Broad_UP",

    "Broad_DOWN",

    "Strict_UP",

    "Strict_DOWN",

    "Extended_signature_total",

    "Core_UP",

    "Core_DOWN",

    "Core_signature_total",

    "UP_retention_percent",

    "DOWN_retention_percent",

    "Immune_genes_flagged"

  ),

  value = c(

    length(broad_up_genes),

    length(broad_down_genes),

    length(strict_up_genes),

    length(strict_down_genes),

    nrow(extended_signature),

    nrow(core_up),

    nrow(core_down),

    nrow(core_signature),

    round(
      up_retention,
      2
    ),

    round(
      down_retention,
      2
    ),

    sum(
      ranked_signature$immune_flag
    )

  )

)


write.csv(
  final_summary,
  file.path(
    output_dir,
    "HER2_Luminal_FINAL_SIGNATURE_SUMMARY.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. FINAL REPORT
# ============================================================

cat("\n")
cat("==============================================\n")
cat("SCRIPT 29 COMPLETE\n")
cat("==============================================\n\n")


cat(
  "Extended UP:",
  nrow(extended_up),
  "\n"
)

cat(
  "Extended DOWN:",
  nrow(extended_down),
  "\n"
)

cat(
  "Extended TOTAL:",
  nrow(extended_signature),
  "\n\n"
)


cat(
  "Core UP:",
  nrow(core_up),
  "\n"
)

cat(
  "Core DOWN:",
  nrow(core_down),
  "\n"
)

cat(
  "Core TOTAL:",
  nrow(core_signature),
  "\n\n"
)


cat(
  "UP retention:",
  round(
    up_retention,
    2
  ),
  "%\n"
)

cat(
  "DOWN retention:",
  round(
    down_retention,
    2
  ),
  "%\n"
)

cat(
  "Immune genes flagged:",
  sum(
    ranked_signature$immune_flag
  ),
  "\n\n"
)


cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)


cat("KEY FINAL FILES:\n\n")

cat(
  "1. HER2_Luminal_CORE_SIGNATURE.csv\n"
)

cat(
  "2. HER2_Luminal_EXTENDED_SIGNATURE.csv\n"
)

cat(
  "3. HER2_Luminal_FINAL_RANKED_SIGNATURE.csv\n"
)

cat(
  "4. contextual_HER2_luminal_marker_QC.csv\n"
)

cat(
  "5. HER2_Luminal_signature_GO_pathway_membership.csv\n"
)

cat(
  "6. HER2_Luminal_FINAL_SIGNATURE_SUMMARY.csv\n"
)

cat(
  "7. HER2_Luminal_final_signature_QC.png\n\n"
)


cat("IMPORTANT:\n")

cat(
  "The CORE signature is the high-confidence subset.\n"
)

cat(
  "The EXTENDED signature preserves the complete strict signature.\n"
)

cat(
  "The FINAL_RANKED_SIGNATURE preserves gene direction and\n"
)

cat(
  "ranking information for downstream drug-repurposing analysis.\n"
)

cat(
  "No additional biological filtering was performed silently.\n\n"
)


cat(
  "Next step: Script 30 — drug-repurposing / connectivity analysis.\n"
)