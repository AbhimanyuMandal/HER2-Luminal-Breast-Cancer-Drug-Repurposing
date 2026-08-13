# ============================================================
# SCRIPT 36 — FINAL PORTFOLIO CANDIDATE SUMMARY
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Purpose:
# Consolidate computational, mechanistic, and external
# validation evidence into recruiter/portfolio-ready outputs.
#
# This script does NOT claim therapeutic efficacy.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
})

cat("==============================================\n")
cat("SCRIPT 36 — FINAL PORTFOLIO CANDIDATE SUMMARY\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

base_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

ranking_file <- file.path(
    base_dir,
    "HER2_Luminal_FINAL_CANDIDATE_RANKING_QC.csv"
)

validation_file <- file.path(
    base_dir,
    "external_validation",
    "HER2_Luminal_TOP20_EXTERNAL_VALIDATION.csv"
)

out_dir <- file.path(
    base_dir,
    "final_portfolio"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. INPUT CHECK
# ------------------------------------------------------------

cat("Checking required files...\n")

if (!file.exists(ranking_file)) {
    stop(
        "Final ranking file not found:\n",
        ranking_file
    )
}

if (!file.exists(validation_file)) {
    stop(
        "External validation file not found:\n",
        validation_file
    )
}

cat("All required files found.\n\n")


# ------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------

ranking <- read_csv(
    ranking_file,
    show_col_types = FALSE
)

validation <- read_csv(
    validation_file,
    show_col_types = FALSE
)

cat(
    "Final ranking candidates:",
    nrow(ranking),
    "\n"
)

cat(
    "Externally validated candidates:",
    nrow(validation),
    "\n\n"
)


# ------------------------------------------------------------
# 4. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_ranking_cols <- c(
    "final_rank",
    "pert_id",
    "pert_iname",
    "robust_rank",
    "tier",
    "cell_lines_tested",
    "positive_cell_lines",
    "negative_cell_lines",
    "median_score",
    "mean_score",
    "min_score",
    "max_score",
    "reproducibility_score",
    "positive_fraction",
    "chembl_id",
    "pref_name",
    "identity_confidence",
    "mechanism_status_corrected",
    "mechanism_record_count",
    "unique_targets",
    "target_chembl_ids",
    "mechanisms",
    "action_types",
    "final_evidence_score",
    "final_priority"
)

missing_cols <- setdiff(
    required_ranking_cols,
    colnames(ranking)
)

if (length(missing_cols) > 0) {

    stop(
        "Missing required ranking columns:\n",
        paste(missing_cols, collapse = ", ")
    )

}


# ------------------------------------------------------------
# 5. MERGE EXTERNAL VALIDATION
# ------------------------------------------------------------

cat("==============================================\n")
cat("MERGING EXTERNAL VALIDATION\n")
cat("==============================================\n\n")

external_cols <- validation %>%

    select(
        pert_id,
        literature_status,
        publication_count,
        oncology_status,
        oncology_publications,
        external_evidence_class,
        overall_evidence_class
    ) %>%

    distinct(
        pert_id,
        .keep_all = TRUE
    )


final_portfolio <- ranking %>%

    left_join(
        external_cols,
        by = "pert_id"
    )


cat(
    "Final portfolio candidates:",
    nrow(final_portfolio),
    "\n"
)

cat(
    "Candidates with external validation:",
    sum(
        !is.na(
            final_portfolio$literature_status
        )
    ),
    "\n\n"
)


# ------------------------------------------------------------
# 6. CREATE PORTFOLIO EVIDENCE CLASS
# ------------------------------------------------------------

final_portfolio <- final_portfolio %>%

    mutate(

        external_validation_status = case_when(

            literature_status ==
                "EVIDENCE_FOUND" ~
                "Literature signal identified",

            literature_status ==
                "NO_MATCH_FOUND" ~
                "No breast-cancer/HER2 literature signal found",

            literature_status ==
                "QUERY_FAILED" ~
                "External validation unavailable",

            TRUE ~
                "Not externally screened"

        ),

        portfolio_evidence_level = case_when(

            final_priority == "PRIORITY_A" &
            mechanism_status_corrected ==
                "annotated" &
            literature_status ==
                "EVIDENCE_FOUND" ~
                "HIGH",

            final_priority %in%
                c("PRIORITY_A", "PRIORITY_B") &
            literature_status ==
                "EVIDENCE_FOUND" ~
                "MODERATE",

            final_priority %in%
                c("PRIORITY_A", "PRIORITY_B") ~
                "COMPUTATIONAL",

            final_priority == "PRIORITY_C" ~
                "SECONDARY",

            TRUE ~
                "EXPLORATORY"

        )

    )


# ------------------------------------------------------------
# 7. CREATE CLEAN PORTFOLIO TABLE
# ------------------------------------------------------------

portfolio_table <- final_portfolio %>%

    select(

        final_rank,
        pert_id,
        pert_iname,

        final_priority,
        portfolio_evidence_level,

        robust_rank,
        tier,

        cell_lines_tested,
        positive_cell_lines,
        negative_cell_lines,

        median_score,
        mean_score,
        min_score,
        max_score,

        reproducibility_score,
        positive_fraction,

        chembl_id,
        pref_name,
        identity_confidence,

        mechanism_status_corrected,
        mechanism_record_count,
        unique_targets,

        target_chembl_ids,
        mechanisms,
        action_types,

        literature_status,
        publication_count,

        oncology_status,
        oncology_publications,

        external_evidence_class,
        overall_evidence_class,

        final_evidence_score

    ) %>%

    arrange(
        final_rank
    )


# ------------------------------------------------------------
# 8. SAVE COMPLETE PORTFOLIO TABLE
# ------------------------------------------------------------

write_csv(
    portfolio_table,
    file.path(
        out_dir,
        "HER2_Luminal_FINAL_PORTFOLIO_CANDIDATES.csv"
    )
)


# ------------------------------------------------------------
# 9. TOP 10
# ------------------------------------------------------------

top10 <- portfolio_table %>%

    slice_head(
        n = 10
    )

write_csv(
    top10,
    file.path(
        out_dir,
        "HER2_Luminal_TOP10_PORTFOLIO_CANDIDATES.csv"
    )
)


# ------------------------------------------------------------
# 10. PROJECT-LEVEL KEY RESULTS
# ------------------------------------------------------------

key_results <- tibble(

    metric = c(

        "Disease signature genes",
        "Core signature genes",
        "LINCS genes mapped",
        "LINCS signatures screened",
        "LINCS compounds represented",
        "Breast cancer cell lines",
        "Compounds tested in all 5 cell lines",
        "Robust LINCS candidates",
        "Mechanism-supported candidates",
        "Final Priority A candidates",
        "Final Priority B candidates",
        "Final Priority C candidates",
        "Exploratory candidates",
        "Top candidates externally screened"

    ),

    value = c(

        454,
        10,
        315,
        9108,
        6720,
        5,
        167,
        66,
        23,
        3,
        8,
        29,
        26,
        20

    ),

    context = c(

        "Extended disease signature",
        "High-confidence core subset",
        "Disease genes successfully mapped to LINCS",
        "24 h / 10 µM compound signatures",
        "Unique LINCS perturbations represented",
        "MCF7, SKBR3, BT20, MDAMB231, HS578T",
        "Compounds with signatures in every selected cell line",
        "Cross-cell-line robust reversal candidates",
        "Candidates with mechanism evidence",
        "Highest-priority computational candidates",
        "Strong secondary candidates",
        "Additional candidates",
        "Lower-confidence/exploratory candidates",
        "Independent literature screening"

    )

)


write_csv(
    key_results,
    file.path(
        out_dir,
        "HER2_Luminal_PROJECT_KEY_RESULTS.csv"
    )
)


# ------------------------------------------------------------
# 11. PRIORITY SUMMARY
# ------------------------------------------------------------

priority_summary <- portfolio_table %>%

    count(
        final_priority,
        name = "candidates"
    ) %>%

    arrange(
        desc(candidates)
    )


write_csv(
    priority_summary,
    file.path(
        out_dir,
        "SCRIPT36_priority_summary.csv"
    )
)


# ------------------------------------------------------------
# 12. EXTERNAL EVIDENCE SUMMARY
# ------------------------------------------------------------

external_summary <- portfolio_table %>%

    filter(
        !is.na(literature_status)
    ) %>%

    count(
        external_evidence_class,
        name = "candidates"
    ) %>%

    arrange(
        desc(candidates)
    )


write_csv(
    external_summary,
    file.path(
        out_dir,
        "SCRIPT36_external_evidence_summary.csv"
    )
)


# ------------------------------------------------------------
# 13. PRINT TOP 10
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP 10 PORTFOLIO CANDIDATES\n")
cat("==============================================\n\n")

print(

    top10 %>%

        select(
            final_rank,
            pert_iname,
            final_priority,
            portfolio_evidence_level,
            positive_cell_lines,
            median_score,
            reproducibility_score,
            identity_confidence,
            mechanism_status_corrected,
            literature_status
        ),

    n = 10
)


# ------------------------------------------------------------
# 14. PRINT KEY PROJECT NUMBERS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("PROJECT KEY RESULTS\n")
cat("==============================================\n\n")

print(
    key_results
)


# ------------------------------------------------------------
# 15. FINAL QC
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("FINAL PORTFOLIO QC\n")
cat("==============================================\n\n")

cat(
    "Final candidates:",
    nrow(portfolio_table),
    "\n"
)

cat(
    "Top 10 candidates:",
    nrow(top10),
    "\n"
)

cat(
    "Priority A:",
    sum(
        portfolio_table$final_priority ==
            "PRIORITY_A",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Priority B:",
    sum(
        portfolio_table$final_priority ==
            "PRIORITY_B",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Priority C:",
    sum(
        portfolio_table$final_priority ==
            "PRIORITY_C",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Exploratory:",
    sum(
        portfolio_table$final_priority ==
            "EXPLORATORY",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Mechanism-supported:",
    sum(
        portfolio_table$mechanism_status_corrected ==
            "annotated",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "External literature screened:",
    sum(
        !is.na(
            portfolio_table$literature_status
        )
    ),
    "\n\n"
)


cat("==============================================\n")
cat("SCRIPT 36 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "FINAL ANALYTICAL PIPELINE COMPLETE.\n\n"
)
