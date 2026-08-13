# ============================================================
# SCRIPT 30B — LINCS EXTRACTION QC / METADATA REPAIR
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# IMPORTANT:
# This script does NOT read the 22-GB GCTX file.
# It uses the matrix already extracted by Script 30.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

cat("==============================================\n")
cat("SCRIPT 30B — LINCS EXTRACTION QC\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

lincs_dir <- "/Volumes/T7/LINCS/GSE92742"

out_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

matrix_file <- file.path(
    out_dir,
    "HER2_Luminal_LINCS_315genes_9108signatures.rds"
)

gene_map_file <- file.path(
    out_dir,
    "HER2_Luminal_gene_to_LINCS_mapping.csv"
)

filtered_sig_file <- file.path(
    out_dir,
    "LINCS_primary_filtered_signatures.csv"
)

gene_order_file <- file.path(
    out_dir,
    "HER2_Luminal_LINCS_gene_order.csv"
)


# ------------------------------------------------------------
# 2. CHECK REQUIRED FILES
# ------------------------------------------------------------

required_files <- c(
    matrix_file,
    gene_map_file,
    filtered_sig_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Required Script 30 output(s) missing:\n",
        paste(missing_files, collapse = "\n")
    )

}


# ------------------------------------------------------------
# 3. LOAD EXTRACTED MATRIX
# ------------------------------------------------------------

cat("Loading previously extracted LINCS matrix...\n\n")

lincs_mat <- readRDS(matrix_file)

cat(
    "Matrix dimensions:",
    nrow(lincs_mat),
    "genes x",
    ncol(lincs_mat),
    "signatures\n\n"
)


# ------------------------------------------------------------
# 4. BASIC MATRIX QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("MATRIX QC\n")
cat("==============================================\n\n")

cat(
    "Genes:",
    nrow(lincs_mat),
    "\n"
)

cat(
    "Signatures:",
    ncol(lincs_mat),
    "\n"
)

cat(
    "Missing values:",
    sum(is.na(lincs_mat)),
    "\n"
)

cat(
    "Finite values:",
    sum(is.finite(lincs_mat)),
    "\n"
)

cat(
    "Minimum value:",
    min(lincs_mat, na.rm = TRUE),
    "\n"
)

cat(
    "Maximum value:",
    max(lincs_mat, na.rm = TRUE),
    "\n\n"
)


# ------------------------------------------------------------
# 5. LOAD GENE MAPPING
# ------------------------------------------------------------

cat("Loading LINCS gene mapping...\n\n")

gene_map <- read_csv(
    gene_map_file,
    show_col_types = FALSE
)

cat(
    "Mapped disease genes:",
    nrow(gene_map),
    "\n\n"
)


# ------------------------------------------------------------
# 6. FIX LINCS GENE ID TYPE
# ------------------------------------------------------------

gene_map <- gene_map %>%
    mutate(
        pr_gene_id = as.character(pr_gene_id)
    )


# ------------------------------------------------------------
# 7. CHECK MATRIX GENE IDs
# ------------------------------------------------------------

matrix_gene_ids <- rownames(lincs_mat)

if (is.null(matrix_gene_ids)) {

    stop(
        "The extracted LINCS matrix has no row names."
    )

}

matrix_gene_ids <- as.character(
    matrix_gene_ids
)

cat(
    "Matrix gene IDs:",
    length(matrix_gene_ids),
    "\n"
)

cat(
    "Mapped gene IDs:",
    nrow(gene_map),
    "\n\n"
)


# ------------------------------------------------------------
# 8. CREATE GENE ORDER TABLE
# ------------------------------------------------------------

gene_order <- data.frame(
    LINCS_gene_id = matrix_gene_ids,
    stringsAsFactors = FALSE
)

gene_order <- gene_order %>%
    left_join(
        gene_map %>%
            select(
                pr_gene_id,
                gene,
                logFC,
                direction
            ),
        by = c(
            "LINCS_gene_id" = "pr_gene_id"
        )
    )


# ------------------------------------------------------------
# 9. GENE MAPPING QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("GENE MAPPING QC\n")
cat("==============================================\n\n")

cat(
    "Genes extracted:",
    nrow(gene_order),
    "\n"
)

cat(
    "Genes with disease-symbol mapping:",
    sum(!is.na(gene_order$gene)),
    "\n"
)

cat(
    "Genes without disease-symbol mapping:",
    sum(is.na(gene_order$gene)),
    "\n\n"
)

cat("Example gene mapping:\n\n")

print(
    head(
        gene_order,
        20
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 10. LOAD FILTERED SIGNATURE METADATA
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SIGNATURE METADATA QC\n")
cat("==============================================\n\n")

filtered_sig <- read_csv(
    filtered_sig_file,
    show_col_types = FALSE
)

cat(
    "Filtered signatures:",
    nrow(filtered_sig),
    "\n"
)

cat(
    "Unique signature IDs:",
    n_distinct(filtered_sig$sig_id),
    "\n"
)

cat(
    "Unique compounds:",
    n_distinct(filtered_sig$pert_id),
    "\n\n"
)


# ------------------------------------------------------------
# 11. CHECK MATRIX SIGNATURE IDS
# ------------------------------------------------------------

matrix_sig_ids <- colnames(lincs_mat)

if (is.null(matrix_sig_ids)) {

    stop(
        "The extracted LINCS matrix has no column names."
    )

}

matrix_sig_ids <- as.character(
    matrix_sig_ids
)

cat(
    "Signature IDs in matrix:",
    length(matrix_sig_ids),
    "\n"
)


# ------------------------------------------------------------
# 12. CHECK SIGNATURE OVERLAP
# ------------------------------------------------------------

signature_overlap <- intersect(
    matrix_sig_ids,
    filtered_sig$sig_id
)

cat(
    "Signature IDs matching metadata:",
    length(signature_overlap),
    "\n"
)

cat(
    "Signature IDs missing from metadata:",
    length(matrix_sig_ids) -
        length(signature_overlap),
    "\n\n"
)


# ------------------------------------------------------------
# 13. CHECK CELL-LINE DISTRIBUTION
# ------------------------------------------------------------

cat("==============================================\n")
cat("CELL-LINE DISTRIBUTION\n")
cat("==============================================\n\n")

print(
    filtered_sig %>%
        count(
            cell_id,
            sort = TRUE
        )
)


# ------------------------------------------------------------
# 14. CHECK TREATMENT CONDITIONS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TREATMENT CONDITIONS\n")
cat("==============================================\n\n")

cat(
    "Treatment times:\n"
)

print(
    table(
        filtered_sig$pert_itime,
        useNA = "ifany"
    )
)

cat("\nDose labels:\n")

print(
    table(
        filtered_sig$pert_idose,
        useNA = "ifany"
    )
)


# ------------------------------------------------------------
# 15. SAVE REPAIRED GENE ORDER
# ------------------------------------------------------------

write_csv(
    gene_order,
    gene_order_file
)


# ------------------------------------------------------------
# 16. SAVE FINAL EXTRACTION SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(

    disease_signature_genes = 454,

    LINCS_mapped_genes = nrow(
        gene_map
    ),

    genes_extracted = nrow(
        lincs_mat
    ),

    signatures_extracted = ncol(
        lincs_mat
    ),

    signatures_in_metadata = nrow(
        filtered_sig
    ),

    signature_metadata_matches = length(
        signature_overlap
    ),

    missing_matrix_values = sum(
        is.na(lincs_mat)
    ),

    finite_matrix_values = sum(
        is.finite(lincs_mat)
    ),

    stringsAsFactors = FALSE
)

write_csv(
    summary_df,
    file.path(
        out_dir,
        "SCRIPT30B_extraction_QC_summary.csv"
    )
)


# ------------------------------------------------------------
# 17. FINAL STATUS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 30B COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Extracted matrix:",
    nrow(lincs_mat),
    "genes x",
    ncol(lincs_mat),
    "signatures\n"
)

cat(
    "Gene metadata matches:",
    sum(!is.na(gene_order$gene)),
    "\n"
)

cat(
    "Signature metadata matches:",
    length(signature_overlap),
    "\n"
)

cat(
    "Missing matrix values:",
    sum(is.na(lincs_mat)),
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "The large GCTX file was NOT accessed.\n"
)

cat(
    "NEXT STEP: Script 31 — LINCS connectivity scoring\n"
)