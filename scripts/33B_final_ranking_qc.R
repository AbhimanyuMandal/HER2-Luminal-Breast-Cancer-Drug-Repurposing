# ============================================================
# SCRIPT 33B — FINAL RANKING QC & CORRECTION
# HER2+ Luminal Breast Cancer Drug Repurposing
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

cat("==============================================\n")
cat("SCRIPT 33B — FINAL RANKING QC & CORRECTION\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

base_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

ranking_file <- file.path(
    base_dir,
    "HER2_Luminal_FINAL_CANDIDATE_RANKING.csv"
)

annotation_file <- file.path(
    base_dir,
    "HER2_Luminal_candidate_mechanism_annotation.csv"
)

evidence_file <- file.path(
    base_dir,
    "SCRIPT32C_candidate_evidence_annotation.csv"
)

summary_file <- file.path(
    base_dir,
    "SCRIPT33_final_ranking_summary.csv"
)


# ------------------------------------------------------------
# 2. CHECK INPUTS
# ------------------------------------------------------------

cat("Checking input files...\n")

required_files <- c(
    ranking_file,
    annotation_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing required files:\n",
        paste(missing_files, collapse = "\n")
    )

}

cat("Required files found.\n\n")


# ------------------------------------------------------------
# 3. LOAD FINAL RANKING
# ------------------------------------------------------------

cat("Loading Script 33 final ranking...\n")

ranking <- read_csv(
    ranking_file,
    show_col_types = FALSE
)

cat(
    "Ranking rows:",
    nrow(ranking),
    "\n"
)

cat(
    "Ranking columns:",
    ncol(ranking),
    "\n\n"
)


# ------------------------------------------------------------
# 4. LOAD ANNOTATION
# ------------------------------------------------------------

cat("Loading candidate annotation...\n")

annotation <- read_csv(
    annotation_file,
    show_col_types = FALSE
)

cat(
    "Annotation rows:",
    nrow(annotation),
    "\n"
)

cat(
    "Annotation columns:",
    ncol(annotation),
    "\n\n"
)


# ------------------------------------------------------------
# 5. INPUT STRUCTURE QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("INPUT STRUCTURE QC\n")
cat("==============================================\n\n")

cat("Ranking columns:\n")
print(colnames(ranking))

cat("\nAnnotation columns:\n")
print(colnames(annotation))


required_ranking_cols <- c(
    "pert_id",
    "pert_iname",
    "robust_rank",
    "tier",
    "positive_cell_lines",
    "median_score",
    "mean_score",
    "reproducibility_score",
    "identity_confidence",
    "mechanism_status_corrected",
    "final_evidence_score",
    "final_priority"
)

missing_ranking_cols <- setdiff(
    required_ranking_cols,
    colnames(ranking)
)

if (length(missing_ranking_cols) > 0) {

    stop(
        "Missing required ranking columns:\n",
        paste(missing_ranking_cols, collapse = ", ")
    )

}

cat("\nRanking structure QC: PASSED\n")


# ------------------------------------------------------------
# 6. RECOVER CELL-LINE INFORMATION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("CELL-LINE QC CORRECTION\n")
cat("==============================================\n\n")

# Script 33 ranking does not contain cell_id.
# Recover cell-line information from Script 31 connectivity results.

connectivity_file <- file.path(
    base_dir,
    "HER2_Luminal_LINCS_cell_line_connectivity.csv"
)

if (file.exists(connectivity_file)) {

    connectivity <- read_csv(
        connectivity_file,
        show_col_types = FALSE
    )

    if ("cell_id" %in% colnames(connectivity)) {

        represented_cell_lines <- sort(
            unique(connectivity$cell_id)
        )

        cat(
            "Cell lines represented:\n",
            paste(
                represented_cell_lines,
                collapse = ", "
            ),
            "\n"
        )

        cat(
            "Number of cell lines:",
            length(represented_cell_lines),
            "\n"
        )

        cell_line_qc_status <- "PASSED"

    } else {

        warning(
            "cell_id not found in connectivity file."
        )

        represented_cell_lines <- character(0)
        cell_line_qc_status <- "NOT_AVAILABLE"

    }

} else {

    warning(
        "Connectivity file not found. Cell-line QC cannot be reconstructed."
    )

    represented_cell_lines <- character(0)
    cell_line_qc_status <- "NOT_AVAILABLE"

}


# ------------------------------------------------------------
# 7. COMPOUND-LEVEL UNIQUENESS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("COMPOUND UNIQUENESS QC\n")
cat("==============================================\n\n")

unique_compounds <- n_distinct(
    ranking$pert_id
)

cat(
    "Unique candidate compounds:",
    unique_compounds,
    "\n"
)

duplicate_pert_ids <- ranking %>%
    count(pert_id) %>%
    filter(n > 1)

cat(
    "Duplicated LINCS compound IDs:",
    nrow(duplicate_pert_ids),
    "\n"
)

if (nrow(duplicate_pert_ids) == 0) {

    cat("Compound uniqueness QC: PASSED\n")

} else {

    cat(
        "WARNING:",
        nrow(duplicate_pert_ids),
        "compound IDs occur more than once.\n"
    )

}


# ------------------------------------------------------------
# 8. CORRECT max_phase VALUES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("ChEMBL PHASE QC\n")
cat("==============================================\n\n")

if ("max_phase" %in% colnames(annotation)) {

    phase_raw <- annotation$max_phase

    phase_numeric <- suppressWarnings(
        as.numeric(phase_raw)
    )

    phase_numeric[
        is.infinite(phase_numeric)
    ] <- NA_real_

    annotation$max_phase_corrected <- phase_numeric

    cat(
        "Missing/invalid max_phase values:",
        sum(is.na(phase_numeric)),
        "\n"
    )

    cat(
        "Valid max_phase values:",
        sum(!is.na(phase_numeric)),
        "\n"
    )

} else {

    annotation$max_phase_corrected <- NA_real_

    cat(
        "max_phase column unavailable.\n"
    )

}


# ------------------------------------------------------------
# 9. JOIN CLEANED ANNOTATION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("ANNOTATION CONSOLIDATION\n")
cat("==============================================\n\n")

annotation_unique <- annotation %>%
    group_by(pert_id) %>%
    summarise(

        chembl_id = first(
            na.omit(chembl_id),
            default = NA_character_
        ),

        pref_name = first(
            na.omit(pref_name),
            default = NA_character_
        ),

        max_phase = if (
            all(is.na(max_phase_corrected))
        ) {
            NA_real_
        } else {
            max(
                max_phase_corrected,
                na.rm = TRUE
            )
        },

        match_status = first(
            na.omit(match_status),
            default = NA_character_
        ),

        target_chembl_id = paste(
            unique(
                na.omit(target_chembl_id)
            ),
            collapse = ";"
        ),

        mechanism_of_action = paste(
            unique(
                na.omit(mechanism_of_action)
            ),
            collapse = ";"
        ),

        action_type = paste(
            unique(
                na.omit(action_type)
            ),
            collapse = ";"
        ),

        .groups = "drop"

    )


cat(
    "Unique annotated compounds:",
    nrow(annotation_unique),
    "\n"
)


# ------------------------------------------------------------
# 10. MERGE WITH FINAL RANKING
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("FINAL EVIDENCE TABLE\n")
cat("==============================================\n\n")

final_table <- ranking %>%

    left_join(
        annotation_unique,
        by = "pert_id",
        suffix = c("", "_annotation")
    ) %>%

    arrange(
        final_rank
    )


# ------------------------------------------------------------
# 11. FINAL MISSING-VALUE QC
# ------------------------------------------------------------

cat("Connectivity missing:",
    sum(is.na(final_table$median_score)),
    "\n"
)

cat(
    "Identity confidence missing:",
    sum(is.na(final_table$identity_confidence)),
    "\n"
)

cat(
    "Mechanism status missing:",
    sum(is.na(
        final_table$mechanism_status_corrected
    )),
    "\n"
)

cat(
    "Final evidence score missing:",
    sum(is.na(
        final_table$final_evidence_score
    )),
    "\n"
)


# ------------------------------------------------------------
# 12. PRIORITY DISTRIBUTION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("FINAL PRIORITY DISTRIBUTION\n")
cat("==============================================\n\n")

priority_distribution <- final_table %>%

    count(
        final_priority,
        sort = FALSE
    )

print(priority_distribution)


# ------------------------------------------------------------
# 13. PRIORITY A
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("PRIORITY A CANDIDATES\n")
cat("==============================================\n\n")

priority_a <- final_table %>%

    filter(
        final_priority == "PRIORITY_A"
    ) %>%

    arrange(
        final_rank
    )

print(
    priority_a %>%
        select(
            final_rank,
            pert_iname,
            pert_id,
            robust_rank,
            positive_cell_lines,
            median_score,
            reproducibility_score,
            identity_confidence,
            mechanism_status_corrected,
            final_evidence_score,
            chembl_id,
            pref_name,
            mechanism_of_action,
            action_type
        ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 14. HIGH-CONFIDENCE IDENTITY CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("HIGH-CONFIDENCE IDENTITY CANDIDATES\n")
cat("==============================================\n\n")

high_identity <- final_table %>%

    filter(
        identity_confidence == "HIGH"
    ) %>%

    arrange(
        final_rank
    )

cat(
    "High-confidence identities:",
    nrow(high_identity),
    "\n\n"
)

print(
    high_identity %>%
        select(
            final_rank,
            pert_iname,
            pert_id,
            chembl_id,
            pref_name,
            median_score,
            positive_cell_lines,
            final_priority
        ) %>%
        head(30),
    row.names = FALSE
)


# ------------------------------------------------------------
# 15. MECHANISM-SUPPORTED CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("MECHANISM-SUPPORTED CANDIDATES\n")
cat("==============================================\n\n")

mechanism_supported <- final_table %>%

    filter(
        mechanism_status_corrected == "annotated"
    ) %>%

    arrange(
        final_rank
    )

cat(
    "Candidates with mechanism evidence:",
    nrow(mechanism_supported),
    "\n\n"
)

print(
    mechanism_supported %>%
        select(
            final_rank,
            pert_iname,
            pert_id,
            chembl_id,
            pref_name,
            median_score,
            positive_cell_lines,
            final_priority
        ) %>%
        head(30),
    row.names = FALSE
)


# ------------------------------------------------------------
# 16. STRONG LINCS / NO MECHANISM
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("STRONG LINCS REVERSALS WITHOUT MECHANISM\n")
cat("==============================================\n\n")

strong_no_mechanism <- final_table %>%

    filter(
        mechanism_status_corrected != "annotated"
    ) %>%

    arrange(
        final_rank
    )

cat(
    "Candidates:",
    nrow(strong_no_mechanism),
    "\n\n"
)

print(
    strong_no_mechanism %>%
        select(
            final_rank,
            pert_iname,
            pert_id,
            median_score,
            positive_cell_lines,
            reproducibility_score,
            identity_confidence,
            final_priority
        ) %>%
        head(20),
    row.names = FALSE
)


# ------------------------------------------------------------
# 17. DUPLICATE ChEMBL QC
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("DUPLICATE ChEMBL ASSIGNMENTS\n")
cat("==============================================\n\n")

duplicate_chembl <- final_table %>%

    filter(
        !is.na(chembl_id),
        chembl_id != ""
    ) %>%

    group_by(
        chembl_id
    ) %>%

    summarise(
        n_compounds = n_distinct(pert_id),
        compounds = paste(
            unique(pert_iname),
            collapse = "; "
        ),
        .groups = "drop"
    ) %>%

    filter(
        n_compounds > 1
    )

print(
    duplicate_chembl,
    row.names = FALSE
)


# ------------------------------------------------------------
# 18. SAVE CLEAN FINAL TABLE
# ------------------------------------------------------------

clean_final_file <- file.path(
    base_dir,
    "HER2_Luminal_FINAL_CANDIDATE_RANKING_QC.csv"
)

write_csv(
    final_table,
    clean_final_file
)


# ------------------------------------------------------------
# 19. SAVE PRIORITY A
# ------------------------------------------------------------

write_csv(
    priority_a,
    file.path(
        base_dir,
        "HER2_Luminal_PRIORITY_A_CANDIDATES_QC.csv"
    )
)


# ------------------------------------------------------------
# 20. SAVE HIGH-CONFIDENCE IDENTITIES
# ------------------------------------------------------------

write_csv(
    high_identity,
    file.path(
        base_dir,
        "HER2_Luminal_HIGH_CONFIDENCE_IDENTITIES.csv"
    )
)


# ------------------------------------------------------------
# 21. SAVE MECHANISM-SUPPORTED CANDIDATES
# ------------------------------------------------------------

write_csv(
    mechanism_supported,
    file.path(
        base_dir,
        "HER2_Luminal_MECHANISM_SUPPORTED_CANDIDATES.csv"
    )
)


# ------------------------------------------------------------
# 22. SAVE STRONG NO-MECHANISM CANDIDATES
# ------------------------------------------------------------

write_csv(
    strong_no_mechanism,
    file.path(
        base_dir,
        "HER2_Luminal_STRONG_LINCS_NO_MECHANISM_QC.csv"
    )
)


# ------------------------------------------------------------
# 23. FINAL SUMMARY
# ------------------------------------------------------------

final_summary <- data.frame(

    total_candidates =
        nrow(final_table),

    unique_compounds =
        n_distinct(final_table$pert_id),

    priority_A =
        sum(
            final_table$final_priority ==
                "PRIORITY_A",
            na.rm = TRUE
        ),

    priority_B =
        sum(
            final_table$final_priority ==
                "PRIORITY_B",
            na.rm = TRUE
        ),

    priority_C =
        sum(
            final_table$final_priority ==
                "PRIORITY_C",
            na.rm = TRUE
        ),

    exploratory =
        sum(
            final_table$final_priority ==
                "EXPLORATORY",
            na.rm = TRUE
        ),

    high_confidence_identities =
        sum(
            final_table$identity_confidence ==
                "HIGH",
            na.rm = TRUE
        ),

    review_identities =
        sum(
            final_table$identity_confidence ==
                "REVIEW",
            na.rm = TRUE
        ),

    mechanism_supported =
        sum(
            final_table$mechanism_status_corrected ==
                "annotated",
            na.rm = TRUE
        ),

    no_mechanism =
        sum(
            final_table$mechanism_status_corrected !=
                "annotated",
            na.rm = TRUE
        ),

    duplicate_chembl_assignments =
        nrow(duplicate_chembl),

    cell_lines_represented =
        length(represented_cell_lines),

    missing_connectivity =
        sum(
            is.na(final_table$median_score)
        ),

    cell_line_qc_status =
        cell_line_qc_status

)

write_csv(
    final_summary,
    file.path(
        base_dir,
        "SCRIPT33B_final_QC_summary.csv"
    )
)


# ------------------------------------------------------------
# 24. FINAL REPORT
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 33B COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Total candidates:",
    final_summary$total_candidates,
    "\n"
)

cat(
    "Priority A:",
    final_summary$priority_A,
    "\n"
)

cat(
    "Priority B:",
    final_summary$priority_B,
    "\n"
)

cat(
    "Priority C:",
    final_summary$priority_C,
    "\n"
)

cat(
    "Exploratory:",
    final_summary$exploratory,
    "\n\n"
)

cat(
    "High-confidence identities:",
    final_summary$high_confidence_identities,
    "\n"
)

cat(
    "Identity review:",
    final_summary$review_identities,
    "\n"
)

cat(
    "Mechanism-supported:",
    final_summary$mechanism_supported,
    "\n"
)

cat(
    "No mechanism record:",
    final_summary$no_mechanism,
    "\n"
)

cat(
    "Duplicate ChEMBL assignments:",
    final_summary$duplicate_chembl_assignments,
    "\n"
)

cat(
    "Cell lines represented:",
    final_summary$cell_lines_represented,
    "\n"
)

cat(
    "Missing connectivity:",
    final_summary$missing_connectivity,
    "\n\n"
)

cat(
    "Output directory:\n",
    base_dir,
    "\n\n"
)

cat("KEY OUTPUTS:\n")
cat("1. HER2_Luminal_FINAL_CANDIDATE_RANKING_QC.csv\n")
cat("2. HER2_Luminal_PRIORITY_A_CANDIDATES_QC.csv\n")
cat("3. HER2_Luminal_HIGH_CONFIDENCE_IDENTITIES.csv\n")
cat("4. HER2_Luminal_MECHANISM_SUPPORTED_CANDIDATES.csv\n")
cat("5. HER2_Luminal_STRONG_LINCS_NO_MECHANISM_QC.csv\n")
cat("6. SCRIPT33B_final_QC_summary.csv\n")

cat("\n==============================================\n")
cat("NO LINCS GCTX RE-EXTRACTION REQUIRED\n")
cat("==============================================\n")