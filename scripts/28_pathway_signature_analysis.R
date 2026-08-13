# ============================================================
# SCRIPT 28 — PATHWAY / SIGNATURE ANALYSIS
# ============================================================
#
# HER2+ LUMINAL EPITHELIAL DISEASE SIGNATURE
#
# Input:
#   Script 27 outputs
#
# Primary biological question:
#   What pathways and biological programs distinguish
#   HER2+ luminal breast cancer from normal luminal breast?
#
# Analyses:
#   1. GO Biological Process enrichment — UP
#   2. GO Biological Process enrichment — DOWN
#   3. KEGG enrichment — UP
#   4. KEGG enrichment — DOWN
#   5. GSEA using complete ranked DE list
#   6. Broad vs strict signature comparison
#
# IMPORTANT:
#   This script does NOT perform differential expression.
#   This script does NOT alter Script 27 outputs.
#
# ============================================================


# ============================================================
# 0. PACKAGE CHECK
# ============================================================

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ggplot2",
  "dplyr"
)

missing_packages <- required_packages[
  !sapply(
    required_packages,
    requireNamespace,
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {

  stop(
    paste0(
      "\nMissing required R packages:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nInstall them before running Script 28.\n"
    )
  )
}


suppressPackageStartupMessages({

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)

})


cat("==============================================\n")
cat("SCRIPT 28 — PATHWAY / SIGNATURE ANALYSIS\n")
cat("==============================================\n\n")


# ============================================================
# 1. FILE PATHS
# ============================================================

input_dir <- "results/disease_signatures/HER2_Luminal"

output_dir <- paste0(
  "results/pathway_analysis/HER2_Luminal"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. INPUT FILES
# ============================================================

strict_up_file <- file.path(
  input_dir,
  "HER2_Luminal_strict_UP.csv"
)

strict_down_file <- file.path(
  input_dir,
  "HER2_Luminal_strict_DOWN.csv"
)

broad_up_file <- file.path(
  input_dir,
  "HER2_Luminal_broad_UP.csv"
)

broad_down_file <- file.path(
  input_dir,
  "HER2_Luminal_broad_DOWN.csv"
)

all_de_file <- file.path(
  input_dir,
  "HER2_Luminal_significant_DE_all.csv"
)

strict_signature_file <- file.path(
  input_dir,
  "HER2_Luminal_strict_signature.csv"
)


input_files <- c(
  strict_up_file,
  strict_down_file,
  broad_up_file,
  broad_down_file,
  all_de_file,
  strict_signature_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {

  stop(
    paste0(
      "\nMissing Script 27 output(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      "\n\nRun Script 27 first.\n"
    )
  )

}


# ============================================================
# 3. LOAD SIGNATURES
# ============================================================

cat("Loading Script 27 outputs...\n\n")

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

strict_signature <- read.csv(
  strict_signature_file,
  stringsAsFactors = FALSE
)


cat(
  "Strict UP genes:",
  nrow(strict_up),
  "\n"
)

cat(
  "Strict DOWN genes:",
  nrow(strict_down),
  "\n"
)

cat(
  "Broad UP genes:",
  nrow(broad_up),
  "\n"
)

cat(
  "Broad DOWN genes:",
  nrow(broad_down),
  "\n"
)

cat(
  "Significant DE genes:",
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
# 5. GENE SYMBOL → ENTREZ ID CONVERSION
# ============================================================

cat("Converting gene symbols to Entrez IDs...\n\n")


convert_to_entrez <- function(
  genes
) {

  if (length(genes) == 0) {
    return(
      data.frame(
        SYMBOL = character(),
        ENTREZID = character(),
        stringsAsFactors = FALSE
      )
    )
  }

  mapped <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    columns = c(
      "SYMBOL",
      "ENTREZID"
    )
  )

  mapped <- mapped[
    !is.na(mapped$ENTREZID),
    ,
    drop = FALSE
  ]

  mapped <- mapped[
    !duplicated(
      mapped$SYMBOL
    ),
    ,
    drop = FALSE
  ]

  mapped

}


strict_up_map <- convert_to_entrez(
  strict_up_genes
)

strict_down_map <- convert_to_entrez(
  strict_down_genes
)

broad_up_map <- convert_to_entrez(
  broad_up_genes
)

broad_down_map <- convert_to_entrez(
  broad_down_genes
)


cat(
  "Strict UP mapped:",
  nrow(strict_up_map),
  "\n"
)

cat(
  "Strict DOWN mapped:",
  nrow(strict_down_map),
  "\n\n"
)


# ============================================================
# 6. SAVE GENE MAPPING
# ============================================================

write.csv(
  strict_up_map,
  file.path(
    output_dir,
    "HER2_Luminal_strict_UP_gene_mapping.csv"
  ),
  row.names = FALSE
)

write.csv(
  strict_down_map,
  file.path(
    output_dir,
    "HER2_Luminal_strict_DOWN_gene_mapping.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. GO BIOLOGICAL PROCESS ENRICHMENT
# ============================================================

cat("==============================================\n")
cat("GO BIOLOGICAL PROCESS ENRICHMENT\n")
cat("==============================================\n\n")


run_GO <- function(
  entrez_ids,
  direction
) {

  if (length(entrez_ids) < 5) {

    cat(
      "Too few genes for GO:",
      direction,
      "\n"
    )

    return(NULL)

  }

  result <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  if (!is.null(result) &&
      nrow(as.data.frame(result)) > 0) {

    result_df <- as.data.frame(
      result
    )

    write.csv(
      result_df,
      file.path(
        output_dir,
        paste0(
          "HER2_Luminal_GO_BP_",
          direction,
          ".csv"
        )
      ),
      row.names = FALSE
    )

  }

  result

}


go_up <- run_GO(
  strict_up_map$ENTREZID,
  "UP"
)

go_down <- run_GO(
  strict_down_map$ENTREZID,
  "DOWN"
)


# ============================================================
# 8. KEGG ENRICHMENT
# ============================================================

cat("\n==============================================\n")
cat("KEGG PATHWAY ENRICHMENT\n")
cat("==============================================\n\n")


run_KEGG <- function(
  entrez_ids,
  direction
) {

  if (length(entrez_ids) < 5) {

    cat(
      "Too few genes for KEGG:",
      direction,
      "\n"
    )

    return(NULL)

  }

  result <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )

  if (!is.null(result) &&
      nrow(as.data.frame(result)) > 0) {

    result_df <- as.data.frame(
      result
    )

    write.csv(
      result_df,
      file.path(
        output_dir,
        paste0(
          "HER2_Luminal_KEGG_",
          direction,
          ".csv"
        )
      ),
      row.names = FALSE
    )

  }

  result

}


kegg_up <- run_KEGG(
  strict_up_map$ENTREZID,
  "UP"
)

kegg_down <- run_KEGG(
  strict_down_map$ENTREZID,
  "DOWN"
)


# ============================================================
# 9. GSEA — COMPLETE RANKED DE LIST
# ============================================================
#
# Unlike over-representation analysis, GSEA uses the entire
# ranked DE list rather than an arbitrary significance cutoff.
#
# Ranking:
#
#   sign(logFC) * -log10(PValue)
#
# This gives strongly significant genes with large effects
# higher ranking.
#
# ============================================================

cat("\n==============================================\n")
cat("GSEA — COMPLETE RANKED DE LIST\n")
cat("==============================================\n\n")


gsea_data <- all_de %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(logFC),
    !is.na(PValue),
    PValue > 0
  ) %>%
  mutate(
    rank_score =
      sign(logFC) *
      (-log10(PValue))
  ) %>%
  group_by(gene) %>%
  slice_max(
    order_by = abs(rank_score),
    n = 1
  ) %>%
  ungroup()


# Map symbols to Entrez IDs

gsea_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gsea_data$gene),
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "ENTREZID"
  )
)

gsea_map <- gsea_map[
  !is.na(gsea_map$ENTREZID),
  ,
  drop = FALSE
]

gsea_map <- gsea_map[
  !duplicated(
    gsea_map$SYMBOL
  ),
  ,
  drop = FALSE
]


gsea_data <- gsea_data %>%
  inner_join(
    gsea_map,
    by = c(
      "gene" = "SYMBOL"
    )
  )


gene_list <- gsea_data$rank_score

names(gene_list) <- gsea_data$ENTREZID


# Remove duplicate Entrez IDs

gene_list <- gene_list[
  !duplicated(
    names(gene_list)
  )
]

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


cat(
  "Genes in GSEA ranking:",
  length(gene_list),
  "\n\n"
)


# ============================================================
# 10. GSEA GO BIOLOGICAL PROCESS
# ============================================================

gsea_GO <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)


if (!is.null(gsea_GO)) {

  gsea_GO_df <- as.data.frame(
    gsea_GO
  )

  write.csv(
    gsea_GO_df,
    file.path(
      output_dir,
      "HER2_Luminal_GSEA_GO_BP.csv"
    ),
    row.names = FALSE
  )

}


# ============================================================
# 11. GSEA KEGG
# ============================================================

gsea_KEGG <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)


if (!is.null(gsea_KEGG)) {

  gsea_KEGG_df <- as.data.frame(
    gsea_KEGG
  )

  write.csv(
    gsea_KEGG_df,
    file.path(
      output_dir,
      "HER2_Luminal_GSEA_KEGG.csv"
    ),
    row.names = FALSE
  )

}


# ============================================================
# 12. GSEA TOP PATHWAYS
# ============================================================

if (!is.null(gsea_GO) &&
    nrow(as.data.frame(gsea_GO)) > 0) {

  gsea_GO_df <- as.data.frame(
    gsea_GO
  )

  gsea_GO_sig <- gsea_GO_df %>%
    filter(
      p.adjust < 0.05
    ) %>%
    arrange(
      p.adjust
    )

  write.csv(
    head(
      gsea_GO_sig,
      30
    ),
    file.path(
      output_dir,
      "HER2_Luminal_GSEA_GO_BP_top30.csv"
    ),
    row.names = FALSE
  )

}


if (!is.null(gsea_KEGG) &&
    nrow(as.data.frame(gsea_KEGG)) > 0) {

  gsea_KEGG_df <- as.data.frame(
    gsea_KEGG
  )

  gsea_KEGG_sig <- gsea_KEGG_df %>%
    filter(
      p.adjust < 0.05
    ) %>%
    arrange(
      p.adjust
    )

  write.csv(
    head(
      gsea_KEGG_sig,
      30
    ),
    file.path(
      output_dir,
      "HER2_Luminal_GSEA_KEGG_top30.csv"
    ),
    row.names = FALSE
  )

}


# ============================================================
# 13. BROAD VS STRICT SIGNATURE OVERLAP
# ============================================================

cat("\n==============================================\n")
cat("BROAD VS STRICT SIGNATURE ROBUSTNESS\n")
cat("==============================================\n\n")


broad_up_set <- unique(
  broad_up$gene
)

strict_up_set <- unique(
  strict_up$gene
)

broad_down_set <- unique(
  broad_down$gene
)

strict_down_set <- unique(
  strict_down$gene
)


up_overlap <- intersect(
  broad_up_set,
  strict_up_set
)

down_overlap <- intersect(
  broad_down_set,
  strict_down_set
)


robustness_table <- data.frame(

  direction = c(
    "UP",
    "DOWN"
  ),

  broad_genes = c(
    length(broad_up_set),
    length(broad_down_set)
  ),

  strict_genes = c(
    length(strict_up_set),
    length(strict_down_set)
  ),

  overlapping_genes = c(
    length(up_overlap),
    length(down_overlap)
  ),

  retained_percentage = c(

    100 *
      length(up_overlap) /
      length(broad_up_set),

    100 *
      length(down_overlap) /
      length(broad_down_set)

  )

)


print(
  robustness_table,
  row.names = FALSE
)


write.csv(
  robustness_table,
  file.path(
    output_dir,
    "HER2_Luminal_broad_vs_strict_robustness.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SAVE OVERLAPPING SIGNATURES
# ============================================================

write.csv(
  data.frame(
    gene = up_overlap
  ),
  file.path(
    output_dir,
    "HER2_Luminal_UP_robust_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = down_overlap
  ),
  file.path(
    output_dir,
    "HER2_Luminal_DOWN_robust_genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. GO DOTPLOTS
# ============================================================

if (!is.null(go_up) &&
    nrow(as.data.frame(go_up)) > 0) {

  p_go_up <- dotplot(
    go_up,
    showCategory = 20
  ) +
    ggtitle(
      "HER2+ Luminal — GO Biological Process UP"
    )

  ggsave(
    file.path(
      output_dir,
      "HER2_Luminal_GO_BP_UP_dotplot.png"
    ),
    p_go_up,
    width = 10,
    height = 8,
    dpi = 300
  )

}


if (!is.null(go_down) &&
    nrow(as.data.frame(go_down)) > 0) {

  p_go_down <- dotplot(
    go_down,
    showCategory = 20
  ) +
    ggtitle(
      "HER2+ Luminal — GO Biological Process DOWN"
    )

  ggsave(
    file.path(
      output_dir,
      "HER2_Luminal_GO_BP_DOWN_dotplot.png"
    ),
    p_go_down,
    width = 10,
    height = 8,
    dpi = 300
  )

}


# ============================================================
# 16. KEGG DOTPLOTS
# ============================================================

if (!is.null(kegg_up) &&
    nrow(as.data.frame(kegg_up)) > 0) {

  p_kegg_up <- dotplot(
    kegg_up,
    showCategory = 15
  ) +
    ggtitle(
      "HER2+ Luminal — KEGG UP"
    )

  ggsave(
    file.path(
      output_dir,
      "HER2_Luminal_KEGG_UP_dotplot.png"
    ),
    p_kegg_up,
    width = 10,
    height = 8,
    dpi = 300
  )

}


if (!is.null(kegg_down) &&
    nrow(as.data.frame(kegg_down)) > 0) {

  p_kegg_down <- dotplot(
    kegg_down,
    showCategory = 15
  ) +
    ggtitle(
      "HER2+ Luminal — KEGG DOWN"
    )

  ggsave(
    file.path(
      output_dir,
      "HER2_Luminal_KEGG_DOWN_dotplot.png"
    ),
    p_kegg_down,
    width = 10,
    height = 8,
    dpi = 300
  )

}


# ============================================================
# 17. SUMMARY
# ============================================================

go_up_n <- 0
go_down_n <- 0
kegg_up_n <- 0
kegg_down_n <- 0
gsea_go_n <- 0
gsea_kegg_n <- 0


if (!is.null(go_up)) {
  go_up_n <- sum(
    as.data.frame(go_up)$p.adjust < 0.05,
    na.rm = TRUE
  )
}

if (!is.null(go_down)) {
  go_down_n <- sum(
    as.data.frame(go_down)$p.adjust < 0.05,
    na.rm = TRUE
  )
}

if (!is.null(kegg_up)) {
  kegg_up_n <- sum(
    as.data.frame(kegg_up)$p.adjust < 0.05,
    na.rm = TRUE
  )
}

if (!is.null(kegg_down)) {
  kegg_down_n <- sum(
    as.data.frame(kegg_down)$p.adjust < 0.05,
    na.rm = TRUE
  )
}

if (!is.null(gsea_GO)) {
  gsea_go_n <- sum(
    as.data.frame(gsea_GO)$p.adjust < 0.05,
    na.rm = TRUE
  )
}

if (!is.null(gsea_KEGG)) {
  gsea_kegg_n <- sum(
    as.data.frame(gsea_KEGG)$p.adjust < 0.05,
    na.rm = TRUE
  )
}


pathway_summary <- data.frame(

  analysis = c(
    "GO_BP_UP",
    "GO_BP_DOWN",
    "KEGG_UP",
    "KEGG_DOWN",
    "GSEA_GO_BP",
    "GSEA_KEGG"
  ),

  significant_pathways_FDR05 = c(
    go_up_n,
    go_down_n,
    kegg_up_n,
    kegg_down_n,
    gsea_go_n,
    gsea_kegg_n
  )

)


write.csv(
  pathway_summary,
  file.path(
    output_dir,
    "HER2_Luminal_pathway_analysis_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. FINAL REPORT
# ============================================================

cat("\n")
cat("==============================================\n")
cat("SCRIPT 28 COMPLETE\n")
cat("==============================================\n\n")

cat(
  "Strict UP genes:",
  length(strict_up_genes),
  "\n"
)

cat(
  "Strict DOWN genes:",
  length(strict_down_genes),
  "\n"
)

cat(
  "GSEA genes:",
  length(gene_list),
  "\n\n"
)

cat(
  "GO BP significant UP pathways:",
  go_up_n,
  "\n"
)

cat(
  "GO BP significant DOWN pathways:",
  go_down_n,
  "\n"
)

cat(
  "KEGG significant UP pathways:",
  kegg_up_n,
  "\n"
)

cat(
  "KEGG significant DOWN pathways:",
  kegg_down_n,
  "\n"
)

cat(
  "GSEA GO BP significant pathways:",
  gsea_go_n,
  "\n"
)

cat(
  "GSEA KEGG significant pathways:",
  gsea_kegg_n,
  "\n\n"
)

cat(
  "UP signature retained:",
  round(
    robustness_table$retained_percentage[
      robustness_table$direction == "UP"
    ],
    2
  ),
  "% after immune exclusion\n"
)

cat(
  "DOWN signature retained:",
  round(
    robustness_table$retained_percentage[
      robustness_table$direction == "DOWN"
    ],
    2
  ),
  "% after immune exclusion\n\n"
)

cat("Output directory:\n")
cat(
  output_dir,
  "\n\n"
)

cat("IMPORTANT:\n")
cat(
  "Pathway enrichment describes biological programs,\n"
)

cat(
  "not individual drug targets.\n"
)

cat(
  "The strongest pathways should be reviewed for\n"
)

cat(
  "biological coherence before drug repurposing.\n"
)

cat("\nNext step: Script 29 — disease-signature validation / drug-repurposing preparation.\n")