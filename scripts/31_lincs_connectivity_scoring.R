# ============================================================
# SCRIPT 31 — LINCS CONNECTIVITY SCORING
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Input:
#   Script 30 extracted LINCS matrix
#   Script 30 gene order
#   Script 30 filtered signature metadata
#
# Output:
#   Signature-level connectivity scores
#   Cell-line-specific compound scores
#   Cross-cell-line compound ranking
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

cat("==============================================\n")
cat("SCRIPT 31 — LINCS CONNECTIVITY SCORING\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

out_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

matrix_file <- file.path(
    out_dir,
    "HER2_Luminal_LINCS_454x_signatures.rds"
)

gene_order_file <- file.path(
    out_dir,
    "HER2_Luminal_LINCS_gene_order.csv"
)

sig_metadata_file <- file.path(
    out_dir,
    "LINCS_primary_filtered_signatures.csv"
)


# ------------------------------------------------------------
# 2. CHECK INPUTS
# ------------------------------------------------------------

required_files <- c(
    matrix_file,
    gene_order_file,
    sig_metadata_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing required Script 30 outputs:\n",
        paste(missing_files, collapse = "\n")
    )

}


# ------------------------------------------------------------
# 3. LOAD MATRIX
# ------------------------------------------------------------

cat("Loading LINCS expression matrix...\n")

lincs_mat <- readRDS(
    matrix_file
)

cat(
    "Matrix:",
    nrow(lincs_mat),
    "genes x",
    ncol(lincs_mat),
    "signatures\n\n"
)


# ------------------------------------------------------------
# 4. LOAD GENE ORDER
# ------------------------------------------------------------

cat("Loading gene metadata...\n")

gene_order <- read_csv(
    gene_order_file,
    show_col_types = FALSE
)

cat(
    "Gene metadata rows:",
    nrow(gene_order),
    "\n\n"
)


# ------------------------------------------------------------
# 5. LOAD SIGNATURE METADATA
# ------------------------------------------------------------

cat("Loading signature metadata...\n")

sig_info <- read_csv(
    sig_metadata_file,
    show_col_types = FALSE
)

cat(
    "Signature metadata rows:",
    nrow(sig_info),
    "\n\n"
)


# ------------------------------------------------------------
# 6. BASIC INPUT QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("INPUT QC\n")
cat("==============================================\n\n")

cat(
    "Matrix genes:",
    nrow(lincs_mat),
    "\n"
)

cat(
    "Matrix signatures:",
    ncol(lincs_mat),
    "\n"
)

cat(
    "Gene metadata:",
    nrow(gene_order),
    "\n"
)

cat(
    "Signature metadata:",
    nrow(sig_info),
    "\n"
)

cat(
    "Missing matrix values:",
    sum(is.na(lincs_mat)),
    "\n\n"
)


# ------------------------------------------------------------
# 7. ALIGN GENE METADATA TO MATRIX
# ------------------------------------------------------------

gene_order$LINCS_gene_id <-
    as.character(gene_order$LINCS_gene_id)

matrix_gene_ids <-
    as.character(rownames(lincs_mat))

gene_index <- match(
    matrix_gene_ids,
    gene_order$LINCS_gene_id
)

if (any(is.na(gene_index))) {

    stop(
        "Some matrix genes could not be matched to gene metadata."
    )

}

gene_order <- gene_order[
    gene_index,
    ,
    drop = FALSE
]


# ------------------------------------------------------------
# 8. CHECK DISEASE DIRECTION
# ------------------------------------------------------------

if (!all(
    gene_order$direction %in%
        c("UP", "DOWN")
)) {

    stop(
        "Gene direction contains unexpected values."
    )

}

up_idx <- which(
    gene_order$direction == "UP"
)

down_idx <- which(
    gene_order$direction == "DOWN"
)

cat(
    "Disease UP genes:",
    length(up_idx),
    "\n"
)

cat(
    "Disease DOWN genes:",
    length(down_idx),
    "\n\n"
)


# ------------------------------------------------------------
# 9. CONNECTIVITY SCORE FUNCTION
# ------------------------------------------------------------
#
# Desired therapeutic reversal:
#
# Disease UP   -> LINCS DOWN
# Disease DOWN -> LINCS UP
#
# Therefore:
#
# score = mean(LINCS expression for DOWN disease genes)
#         -
#         mean(LINCS expression for UP disease genes)
#
# Positive score:
#     perturbation opposes disease signature
#
# Negative score:
#     perturbation reinforces disease signature
#
# This is calculated after standardizing each signature
# across the disease genes, preventing raw expression
# magnitude from dominating.
# ------------------------------------------------------------

calculate_connectivity <- function(
    expression_vector,
    up_idx,
    down_idx
) {

    z <- as.numeric(
        scale(expression_vector)
    )

    up_mean <- mean(
        z[up_idx],
        na.rm = TRUE
    )

    down_mean <- mean(
        z[down_idx],
        na.rm = TRUE
    )

    down_mean - up_mean
}


# ------------------------------------------------------------
# 10. CALCULATE SIGNATURE-LEVEL SCORES
# ------------------------------------------------------------

cat("==============================================\n")
cat("CALCULATING SIGNATURE CONNECTIVITY\n")
cat("==============================================\n\n")

n_signatures <- ncol(lincs_mat)

connectivity_scores <- numeric(
    n_signatures
)

cat(
    "Signatures:",
    n_signatures,
    "\n\n"
)

for (i in seq_len(n_signatures)) {

    connectivity_scores[i] <-
        calculate_connectivity(
            lincs_mat[, i],
            up_idx,
            down_idx
        )

    if (
        i %% 1000 == 0 ||
        i == n_signatures
    ) {

        cat(
            "Processed:",
            i,
            "/",
            n_signatures,
            "\n"
        )

    }

}


# ------------------------------------------------------------
# 11. CREATE SIGNATURE-LEVEL RESULT
# ------------------------------------------------------------

signature_results <- sig_info %>%
    mutate(
        connectivity_score =
            connectivity_scores
    ) %>%
    arrange(
        desc(connectivity_score)
    )


# ------------------------------------------------------------
# 12. ADD SIGNATURE RANK
# ------------------------------------------------------------

signature_results <- signature_results %>%
    mutate(
        connectivity_rank =
            row_number()
    )


# ------------------------------------------------------------
# 13. SAVE SIGNATURE RESULTS
# ------------------------------------------------------------

write_csv(
    signature_results,
    file.path(
        out_dir,
        "HER2_Luminal_LINCS_signature_connectivity.csv"
    )
)


# ------------------------------------------------------------
# 14. CELL-LINE-SPECIFIC RANKING
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("CELL-LINE CONNECTIVITY\n")
cat("==============================================\n\n")

cell_line_results <- signature_results %>%
    group_by(
        cell_id
    ) %>%
    arrange(
        desc(connectivity_score),
        .by_group = TRUE
    ) %>%
    mutate(
        cell_line_rank =
            row_number(),
        cell_line_percentile =
            1 -
            (
                (cell_line_rank - 1) /
                n()
            )
    ) %>%
    ungroup()


# ------------------------------------------------------------
# 15. SAVE CELL-LINE RESULTS
# ------------------------------------------------------------

write_csv(
    cell_line_results,
    file.path(
        out_dir,
        "HER2_Luminal_LINCS_cell_line_connectivity.csv"
    )
)


# ------------------------------------------------------------
# 16. COMPOUND-LEVEL AGGREGATION
# ------------------------------------------------------------
#
# IMPORTANT:
# We aggregate within cell line first.
#
# This prevents MCF7 from dominating simply because
# it has many more signatures.
# ------------------------------------------------------------

compound_cellline <- cell_line_results %>%

    group_by(
        pert_id,
        pert_iname,
        cell_id
    ) %>%

    summarise(

        signatures =
            n(),

        median_connectivity =
            median(
                connectivity_score,
                na.rm = TRUE
            ),

        mean_connectivity =
            mean(
                connectivity_score,
                na.rm = TRUE
            ),

        best_connectivity =
            max(
                connectivity_score,
                na.rm = TRUE
            ),

        worst_connectivity =
            min(
                connectivity_score,
                na.rm = TRUE
            ),

        .groups = "drop"
    )


# ------------------------------------------------------------
# 17. CROSS-CELL-LINE COMPOUND AGGREGATION
# ------------------------------------------------------------

compound_results <- compound_cellline %>%

    group_by(
        pert_id,
        pert_iname
    ) %>%

    summarise(

        cell_lines_tested =
            n_distinct(cell_id),

        median_cellline_score =
            median(
                median_connectivity,
                na.rm = TRUE
            ),

        mean_cellline_score =
            mean(
                median_connectivity,
                na.rm = TRUE
            ),

        best_cellline_score =
            max(
                median_connectivity,
                na.rm = TRUE
            ),

        worst_cellline_score =
            min(
                median_connectivity,
                na.rm = TRUE
            ),

        .groups = "drop"
    ) %>%

    arrange(
        desc(median_cellline_score)
    ) %>%

    mutate(
        compound_rank =
            row_number()
    )


# ------------------------------------------------------------
# 18. SAVE COMPOUND RESULTS
# ------------------------------------------------------------

write_csv(
    compound_cellline,
    file.path(
        out_dir,
        "HER2_Luminal_compound_cellline_connectivity.csv"
    )
)

write_csv(
    compound_results,
    file.path(
        out_dir,
        "HER2_Luminal_compound_connectivity_ranked.csv"
    )
)


# ------------------------------------------------------------
# 19. TOP COMPOUNDS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP COMPOUNDS\n")
cat("==============================================\n\n")

print(
    compound_results %>%
        select(
            compound_rank,
            pert_id,
            pert_iname,
            cell_lines_tested,
            median_cellline_score,
            mean_cellline_score,
            best_cellline_score,
            worst_cellline_score
        ) %>%
        head(30),
    row.names = FALSE
)


# ------------------------------------------------------------
# 20. NEGATIVE / REINFORCING COMPOUNDS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("MOST NEGATIVE COMPOUNDS\n")
cat("==============================================\n\n")

print(
    compound_results %>%
        arrange(
            median_cellline_score
        ) %>%
        select(
            pert_id,
            pert_iname,
            cell_lines_tested,
            median_cellline_score,
            mean_cellline_score,
            best_cellline_score,
            worst_cellline_score
        ) %>%
        head(20),
    row.names = FALSE
)


# ------------------------------------------------------------
# 21. CROSS-CELL-LINE CONSISTENCY
# ------------------------------------------------------------

consistency_results <- compound_cellline %>%

    group_by(
        pert_id,
        pert_iname
    ) %>%

    summarise(

        cell_lines_tested =
            n_distinct(cell_id),

        positive_cell_lines =
            sum(
                median_connectivity > 0
            ),

        negative_cell_lines =
            sum(
                median_connectivity < 0
            ),

        consistent_reversal =
            positive_cell_lines ==
            cell_lines_tested,

        .groups = "drop"
    )


compound_results <- compound_results %>%

    left_join(
        consistency_results,
        by = c(
            "pert_id",
            "pert_iname"
        )
    ) %>%

    arrange(
        desc(
            consistent_reversal
        ),
        desc(
            median_cellline_score
        )
    ) %>%

    mutate(
        robust_rank =
            row_number()
    )


# ------------------------------------------------------------
# 22. SAVE FINAL ROBUST RANKING
# ------------------------------------------------------------

write_csv(
    compound_results,
    file.path(
        out_dir,
        "HER2_Luminal_LINCS_ROBUST_COMPOUND_RANKING.csv"
    )
)


# ------------------------------------------------------------
# 23. SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(

    disease_genes =
        nrow(lincs_mat),

    disease_up_genes =
        length(up_idx),

    disease_down_genes =
        length(down_idx),

    LINCS_signatures =
        ncol(lincs_mat),

    unique_compounds =
        n_distinct(
            sig_info$pert_id
        ),

    cell_lines =
        n_distinct(
            sig_info$cell_id
        ),

    compounds_tested_all_cell_lines =
        sum(
            compound_results$cell_lines_tested == 5
        ),

    stringsAsFactors = FALSE
)

write_csv(
    summary_df,
    file.path(
        out_dir,
        "SCRIPT31_connectivity_summary.csv"
    )
)


# ------------------------------------------------------------
# 24. FINAL STATUS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 31 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Disease genes:",
    nrow(lincs_mat),
    "\n"
)

cat(
    "Disease UP genes:",
    length(up_idx),
    "\n"
)

cat(
    "Disease DOWN genes:",
    length(down_idx),
    "\n"
)

cat(
    "LINCS signatures:",
    ncol(lincs_mat),
    "\n"
)

cat(
    "Unique compounds:",
    n_distinct(sig_info$pert_id),
    "\n"
)

cat(
    "Cell lines:",
    n_distinct(sig_info$cell_id),
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "NEXT STEP:\n",
    "Review connectivity rankings before Script 32.\n"
)