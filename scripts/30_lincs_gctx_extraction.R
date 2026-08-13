# ============================================================
# SCRIPT 30 — LINCS GCTX SIGNATURE EXTRACTION
# HER2+ Luminal Breast Cancer Drug Repurposing
# ============================================================
suppressPackageStartupMessages({
    library(cmapR)
    library(dplyr)
    library(readr)
})

cat("==============================================\n")
cat("SCRIPT 30 — LINCS GCTX EXTRACTION\n")
cat("==============================================\n\n")

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

lincs_dir <- "/Volumes/T7/LINCS/GSE92742"

gctx_file <- file.path(
    lincs_dir,
    "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
)

gene_info_file <- file.path(
    lincs_dir,
    "GSE92742_Broad_LINCS_gene_info.txt.gz"
)

sig_info_file <- file.path(
    lincs_dir,
    "GSE92742_Broad_LINCS_sig_info.txt.gz"
)

out_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. CHECK FILES
# ------------------------------------------------------------

required_files <- c(
    gctx_file,
    gene_info_file,
    sig_info_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing required LINCS files:\n",
        paste(missing_files, collapse = "\n")
    )

}

cat("GCTX file found.\n")
cat("Gene metadata found.\n")
cat("Signature metadata found.\n\n")


# ------------------------------------------------------------
# 3. LOAD LINCS METADATA
# ------------------------------------------------------------

cat("Loading LINCS metadata...\n")

gene_info <- read_tsv(
    gene_info_file,
    show_col_types = FALSE
)

sig_info <- read_tsv(
    sig_info_file,
    show_col_types = FALSE
)

cat(
    "Genes in LINCS:",
    nrow(gene_info),
    "\n"
)

cat(
    "Signatures in LINCS:",
    nrow(sig_info),
    "\n\n"
)


# ------------------------------------------------------------
# 4. LOAD FINAL DISEASE SIGNATURE
# ------------------------------------------------------------

signature_file <- paste0(
    "results/final_disease_signature/",
    "HER2_Luminal/",
    "HER2_Luminal_FINAL_RANKED_SIGNATURE.csv"
)

if (!file.exists(signature_file)) {

    stop(
        "Final ranked disease signature not found:\n",
        signature_file
    )

}

disease_sig <- read_csv(
    signature_file,
    show_col_types = FALSE
)

cat("Disease signature loaded.\n")

cat(
    "Genes:",
    nrow(disease_sig),
    "\n\n"
)


# ------------------------------------------------------------
# 5. VERIFY DISEASE SIGNATURE
# ------------------------------------------------------------

required_sig_cols <- c(
    "gene",
    "logFC",
    "direction"
)

missing_sig_cols <- setdiff(
    required_sig_cols,
    colnames(disease_sig)
)

if (length(missing_sig_cols) > 0) {

    stop(
        "Disease signature is missing columns:\n",
        paste(missing_sig_cols, collapse = ", ")
    )

}


# ------------------------------------------------------------
# 6. RECREATE PRIMARY LINCS FILTER
# ------------------------------------------------------------

breast_lines <- c(
    "MCF7",
    "SKBR3",
    "BT20",
    "MDAMB231",
    "HS578T"
)

compound_sig <- sig_info %>%
    filter(
        pert_type == "trt_cp"
    )

breast_24h_10uM <- compound_sig %>%
    filter(
        cell_id %in% breast_lines,
        pert_itime == "24 h",
        pert_idose == "10 µM"
    )

cat("Primary LINCS filter:\n")
cat("Cell lines:", paste(breast_lines, collapse = ", "), "\n")
cat("Time: 24 h\n")
cat("Dose: 10 µM\n\n")

cat(
    "Selected signatures:",
    nrow(breast_24h_10uM),
    "\n"
)

cat(
    "Unique compounds:",
    n_distinct(breast_24h_10uM$pert_id),
    "\n\n"
)


# ------------------------------------------------------------
# 7. MAP DISEASE GENES TO LINCS GENE IDs
# ------------------------------------------------------------

cat("Mapping disease genes to LINCS gene IDs...\n")

gene_map <- disease_sig %>%
    inner_join(
        gene_info %>%
            select(
                pr_gene_id,
                pr_gene_symbol
            ),
        by = c(
            "gene" = "pr_gene_symbol"
        )
    ) %>%
    mutate(
        pr_gene_id = as.character(pr_gene_id)
    ) %>%
    distinct(
        gene,
        .keep_all = TRUE
    )

cat(
    "Genes mapped to LINCS:",
    nrow(gene_map),
    "\n"
)

cat(
    "Genes not mapped:",
    nrow(disease_sig) - nrow(gene_map),
    "\n\n"
)

if (nrow(gene_map) < 400) {

    warning(
        "Fewer than 400 disease-signature genes mapped to LINCS."
    )

}


# ------------------------------------------------------------
# 8. PREPARE GCTX ROW AND COLUMN IDS
# ------------------------------------------------------------

rid <- as.character(
    gene_map$pr_gene_id
)

cid <- unique(
    breast_24h_10uM$sig_id
)

cat(
    "GCTX genes requested:",
    length(rid),
    "\n"
)

cat(
    "GCTX signatures requested:",
    length(cid),
    "\n\n"
)


# ------------------------------------------------------------
# 9. SAVE FILTER METADATA BEFORE LARGE EXTRACTION
# ------------------------------------------------------------

write_csv(
    breast_24h_10uM,
    file.path(
        out_dir,
        "LINCS_primary_filtered_signatures.csv"
    )
)

write_csv(
    gene_map,
    file.path(
        out_dir,
        "HER2_Luminal_gene_to_LINCS_mapping.csv"
    )
)


# ------------------------------------------------------------
# 10. EXTRACT ONLY REQUIRED GCTX SUBMATRIX
# ------------------------------------------------------------

cat("==============================================\n")
cat("EXTRACTING GCTX SUBMATRIX\n")
cat("==============================================\n\n")

cat("Rows:", length(rid), "\n")
cat("Columns:", length(cid), "\n\n")

gctx_subset <- cmapR::parse_gctx(
    fname = gctx_file,
    rid = rid,
    cid = cid,
    matrix_only = TRUE
)

cat("GCTX extraction complete.\n\n")

lincs_mat <- gctx_subset@mat

cat(
    "Extracted matrix:",
    nrow(lincs_mat),
    "genes x",
    ncol(lincs_mat),
    "signatures\n"
)

# ------------------------------------------------------------
# 11. SAVE MATRIX
# ------------------------------------------------------------

saveRDS(
    lincs_mat,
    file.path(
        out_dir,
        "HER2_Luminal_LINCS_315genes_9108signatures.rds"
    )
)

cat(
    "Saved LINCS expression matrix.\n"
)


# ------------------------------------------------------------
# 12. SAVE GENE ORDER
# ------------------------------------------------------------

gene_order <- data.frame(
    LINCS_gene_id = rownames(lincs_mat),
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

write_csv(
    gene_order,
    file.path(
        out_dir,
        "HER2_Luminal_LINCS_gene_order.csv"
    )
)


# ------------------------------------------------------------
# 13. FINAL QC
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("EXTRACTION QC\n")
cat("==============================================\n\n")

cat(
    "Disease genes requested:",
    length(rid),
    "\n"
)

cat(
    "Genes extracted:",
    nrow(lincs_mat),
    "\n"
)

cat(
    "Signatures requested:",
    length(cid),
    "\n"
)

cat(
    "Signatures extracted:",
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
    "\n\n"
)


# ------------------------------------------------------------
# 14. SAVE SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(
    disease_genes_requested = length(rid),
    disease_genes_extracted = nrow(lincs_mat),
    signatures_requested = length(cid),
    signatures_extracted = ncol(lincs_mat),
    missing_values = sum(is.na(lincs_mat)),
    breast_cell_lines = paste(
        breast_lines,
        collapse = ";"
    ),
    treatment_time = "24 h",
    dose = "10 µM"
)

write_csv(
    summary_df,
    file.path(
        out_dir,
        "SCRIPT30_extraction_summary.csv"
    )
)

cat("==============================================\n")
cat("SCRIPT 30 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "NEXT STEP:\n",
    "Script 31 — connectivity scoring\n"
)