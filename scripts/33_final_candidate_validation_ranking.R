# ============================================================
# SCRIPT 33 — FINAL CANDIDATE VALIDATION & RANKING
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# FINAL INTEGRATION:
#
#   Disease signature
#          ↓
#   LINCS transcriptional reversal
#          ↓
#   Cross-cell-line reproducibility
#          ↓
#   ChEMBL identity confidence
#          ↓
#   Mechanism / target evidence
#          ↓
#   Final candidate prioritization
#
# IMPORTANT:
# This ranking is computational prioritization.
# It is NOT evidence of clinical efficacy.
# ============================================================

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

cat("==============================================\n")
cat("SCRIPT 33 — FINAL CANDIDATE VALIDATION & RANKING\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

base_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

candidate_file <- file.path(
    base_dir,
    "SCRIPT32C_candidate_evidence_annotation.csv"
)

evidence_file <- file.path(
    base_dir,
    "SCRIPT32C_compound_evidence_QC.csv"
)

summary_file <- file.path(
    base_dir,
    "SCRIPT32C_validation_summary.csv"
)

out_dir <- base_dir


# ------------------------------------------------------------
# 2. CHECK INPUTS
# ------------------------------------------------------------

required_files <- c(
    candidate_file,
    evidence_file,
    summary_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing required Script 32C outputs:\n",
        paste(
            missing_files,
            collapse = "\n"
        )
    )

}


# ------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------

cat("Loading Script 32C evidence...\n\n")

candidates <- read_csv(
    candidate_file,
    show_col_types = FALSE
)

evidence <- read_csv(
    evidence_file,
    show_col_types = FALSE
)

cat(
    "Candidate annotation rows:",
    nrow(candidates),
    "\n"
)

cat(
    "Evidence rows:",
    nrow(evidence),
    "\n\n"
)


# ------------------------------------------------------------
# 4. COLLAPSE CANDIDATES TO ONE ROW PER COMPOUND
# ------------------------------------------------------------

cat("==============================================\n")
cat("COMPOUND-LEVEL CONSOLIDATION\n")
cat("==============================================\n\n")

# Candidate annotation may contain multiple rows per compound
# because of multiple mechanism-target records.
#
# The final ranking MUST operate at the compound level.

candidate_compounds <- candidates %>%

    group_by(
        pert_id
    ) %>%

    summarise(

        pert_iname =
            first(na.omit(pert_iname)),

        cell_lines_tested =
            max(cell_lines_tested, na.rm = TRUE),

        positive_cell_lines =
            max(positive_cell_lines, na.rm = TRUE),

        negative_cell_lines =
            max(negative_cell_lines, na.rm = TRUE),

        median_score =
            max(median_score, na.rm = TRUE),

        mean_score =
            max(mean_score, na.rm = TRUE),

        min_score =
            max(min_score, na.rm = TRUE),

        max_score =
            max(max_score, na.rm = TRUE),

        sd_score =
            max(sd_score, na.rm = TRUE),

        positive_fraction =
            max(positive_fraction, na.rm = TRUE),

        negative_fraction =
            max(negative_fraction, na.rm = TRUE),

        reproducibility_score =
            max(reproducibility_score, na.rm = TRUE),

        tier =
            first(na.omit(tier)),

        robust_rank =
            min(robust_rank, na.rm = TRUE),

        canonical_smiles =
            first(
                na.omit(canonical_smiles)
            ),

        inchi_key =
            first(
                na.omit(inchi_key)
            ),

        chembl_id =
            first(
                na.omit(chembl_id)
            ),

        pref_name =
            first(
                na.omit(pref_name)
            ),

        max_phase =
            max(
                max_phase,
                na.rm = TRUE
            ),

        match_status =
            first(
                na.omit(match_status)
            ),

        identity_confidence =
            first(
                na.omit(identity_confidence)
            ),

        mechanism_status_corrected =
            first(
                na.omit(
                    mechanism_status_corrected
                )
            ),

        mechanism_record_count =
            max(
                mechanism_record_count,
                na.rm = TRUE
            ),

        unique_targets =
            max(
                unique_targets,
                na.rm = TRUE
            ),

        target_chembl_ids =
            first(
                na.omit(target_chembl_ids)
            ),

        mechanisms =
            first(
                na.omit(mechanisms)
            ),

        action_types =
            first(
                na.omit(action_types)
            ),

        .groups = "drop"
    )


# Replace invalid Inf values created by max(..., na.rm=TRUE)

candidate_compounds <- candidate_compounds %>%

    mutate(

        across(
            where(is.numeric),
            ~ ifelse(
                is.infinite(.x),
                NA,
                .x
            )
        )

    )


cat(
    "Unique candidate compounds:",
    nrow(candidate_compounds),
    "\n\n"
)


# ------------------------------------------------------------
# 5. MERGE CLEAN EVIDENCE
# ------------------------------------------------------------

candidate_compounds <- candidate_compounds %>%

    left_join(

        evidence %>%

            select(
                pert_id,
                chembl_id,
                identity_confidence,
                mechanism_status_corrected,
                mechanism_record_count,
                unique_targets,
                target_chembl_ids,
                mechanisms,
                action_types
            ) %>%

            distinct(),

        by = c(
            "pert_id",
            "chembl_id"
        ),

        suffix = c(
            "",
            "_evidence"
        )

    ) %>%

    mutate(

        identity_confidence =
            coalesce(
                identity_confidence_evidence,
                identity_confidence
            ),

        mechanism_status_corrected =
            coalesce(
                mechanism_status_corrected_evidence,
                mechanism_status_corrected
            ),

        mechanism_record_count =
            coalesce(
                mechanism_record_count_evidence,
                mechanism_record_count
            ),

        unique_targets =
            coalesce(
                unique_targets_evidence,
                unique_targets
            )

    ) %>%

    select(
        -ends_with("_evidence")
    )


# ------------------------------------------------------------
# 6. INPUT QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("FINAL INPUT QC\n")
cat("==============================================\n\n")

cat(
    "Unique compounds:",
    nrow(candidate_compounds),
    "\n"
)

cat(
    "Cell lines represented:",
    paste(
        sort(
            unique(
                candidates$cell_id[
                    !is.na(candidates$cell_id)
                ]
            )
        ),
        collapse = ", "
    ),
    "\n\n"
)

cat(
    "Missing connectivity scores:",
    sum(
        is.na(
            candidate_compounds$median_score
        )
    ),
    "\n"
)

cat(
    "Missing identity confidence:",
    sum(
        is.na(
            candidate_compounds$identity_confidence
        )
    ),
    "\n"
)

cat(
    "Missing mechanism status:",
    sum(
        is.na(
            candidate_compounds$
                mechanism_status_corrected
        )
    ),
    "\n\n"
)


# ------------------------------------------------------------
# 7. CONNECTIVITY COMPONENT
# ------------------------------------------------------------
#
# Connectivity is the PRIMARY evidence stream.
#
# Normalize using rank among the 66 candidates.
# Higher connectivity = better reversal.
# ------------------------------------------------------------

candidate_compounds <- candidate_compounds %>%

    mutate(

        connectivity_score_norm =
            percent_rank(
                median_score
            ),

        reproducibility_norm =
            percent_rank(
                reproducibility_score
            ),

        positive_fraction_norm =
            percent_rank(
                positive_fraction
            )

    )


# ------------------------------------------------------------
# 8. IDENTITY EVIDENCE
# ------------------------------------------------------------

candidate_compounds <- candidate_compounds %>%

    mutate(

        identity_score =
            case_when(

                identity_confidence == "HIGH" ~
                    1.0,

                identity_confidence == "REVIEW" ~
                    0.5,

                TRUE ~
                    0.0
            )

    )


# ------------------------------------------------------------
# 9. MECHANISM EVIDENCE
# ------------------------------------------------------------
#
# Mechanism evidence is supportive, NOT mandatory.
#
# We do not penalize a compound to zero merely because
# ChEMBL has no mechanism record.
# ------------------------------------------------------------

candidate_compounds <- candidate_compounds %>%

    mutate(

        mechanism_score =
            case_when(

                mechanism_status_corrected ==
                    "annotated" ~ 1.0,

                mechanism_status_corrected ==
                    "no_mechanism_record" ~ 0.0,

                TRUE ~ 0.0
            )

    )


# ------------------------------------------------------------
# 10. FINAL EVIDENCE SCORE
# ------------------------------------------------------------
#
# Weighting:
#
# 50% — transcriptional reversal strength
# 25% — cross-cell-line reproducibility
# 10% — positive-cell-line fraction
# 10% — compound identity confidence
#  5% — mechanism evidence
#
# This prevents ChEMBL availability from overpowering
# the actual LINCS disease-reversal evidence.
# ------------------------------------------------------------

candidate_compounds <- candidate_compounds %>%

    mutate(

        final_evidence_score =

            0.50 *
            connectivity_score_norm +

            0.25 *
            reproducibility_norm +

            0.10 *
            positive_fraction_norm +

            0.10 *
            identity_score +

            0.05 *
            mechanism_score

    )


# ------------------------------------------------------------
# 11. FINAL RANK
# ------------------------------------------------------------

final_ranked <- candidate_compounds %>%

    arrange(
        desc(final_evidence_score),
        desc(median_score),
        desc(reproducibility_score),
        robust_rank
    ) %>%

    mutate(

        final_rank =
            row_number()

    ) %>%

    select(
        final_rank,
        pert_id,
        pert_iname,
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
        max_phase,

        identity_confidence,
        mechanism_status_corrected,
        mechanism_record_count,
        unique_targets,

        target_chembl_ids,
        mechanisms,
        action_types,

        connectivity_score_norm,
        reproducibility_norm,
        positive_fraction_norm,
        identity_score,
        mechanism_score,

        final_evidence_score
    )


# ------------------------------------------------------------
# 12. FINAL TIERS
# ------------------------------------------------------------

final_ranked <- final_ranked %>%

    mutate(

        final_priority = case_when(

            final_rank <= 10 &
            identity_confidence == "HIGH" &
            mechanism_status_corrected ==
                "annotated" ~
                "PRIORITY_A",

            final_rank <= 20 &
            identity_confidence == "HIGH" ~
                "PRIORITY_B",

            final_rank <= 40 ~
                "PRIORITY_C",

            TRUE ~
                "EXPLORATORY"

        )

    )


# ------------------------------------------------------------
# 13. TOP 30
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP 30 FINAL CANDIDATES\n")
cat("==============================================\n\n")

print(

    final_ranked %>%

        select(
            final_rank,
            pert_iname,
            robust_rank,
            tier,
            positive_cell_lines,
            median_score,
            reproducibility_score,
            identity_confidence,
            mechanism_status_corrected,
            final_evidence_score,
            final_priority
        ) %>%

        head(30),

    row.names = FALSE
)


# ------------------------------------------------------------
# 14. PRIORITY DISTRIBUTION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("FINAL PRIORITY DISTRIBUTION\n")
cat("==============================================\n\n")

priority_table <- final_ranked %>%

    count(
        final_priority,
        sort = FALSE
    )

print(
    priority_table
)


# ------------------------------------------------------------
# 15. HIGH-CONFIDENCE FINAL CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("PRIORITY A CANDIDATES\n")
cat("==============================================\n\n")

priority_A <- final_ranked %>%

    filter(
        final_priority == "PRIORITY_A"
    )

print(
    priority_A %>%

        select(
            final_rank,
            pert_iname,
            pert_id,
            median_score,
            positive_cell_lines,
            reproducibility_score,
            chembl_id,
            pref_name,
            mechanisms,
            action_types,
            final_evidence_score
        ),

    row.names = FALSE
)


# ------------------------------------------------------------
# 16. STRONG LINCS REVERSALS WITHOUT MECHANISM
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("STRONG LINCS REVERSALS WITHOUT ChEMBL MECHANISM\n")
cat("==============================================\n\n")

strong_no_mechanism <- final_ranked %>%

    filter(
        mechanism_status_corrected ==
            "no_mechanism_record"
    ) %>%

    arrange(
        final_rank
    ) %>%

    select(
        final_rank,
        pert_iname,
        pert_id,
        median_score,
        positive_cell_lines,
        reproducibility_score,
        identity_confidence,
        chembl_id,
        final_evidence_score
    ) %>%

    head(20)

print(
    strong_no_mechanism,
    row.names = FALSE
)


# ------------------------------------------------------------
# 17. IDENTITY REVIEW CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("HIGH-RANKED IDENTITY-REVIEW CANDIDATES\n")
cat("==============================================\n\n")

identity_review <- final_ranked %>%

    filter(
        identity_confidence == "REVIEW"
    ) %>%

    arrange(
        final_rank
    ) %>%

    select(
        final_rank,
        pert_iname,
        pert_id,
        chembl_id,
        pref_name,
        median_score,
        reproducibility_score,
        final_evidence_score
    ) %>%

    head(20)

print(
    identity_review,
    row.names = FALSE
)


# ------------------------------------------------------------
# 18. DUPLICATE ChEMBL ASSIGNMENTS
# ------------------------------------------------------------

duplicate_final <- final_ranked %>%

    filter(
        !is.na(chembl_id)
    ) %>%

    group_by(
        chembl_id
    ) %>%

    filter(
        n() > 1
    ) %>%

    ungroup()

cat("\n==============================================\n")
cat("DUPLICATE ChEMBL ASSIGNMENTS IN FINAL SET\n")
cat("==============================================\n\n")

if (nrow(duplicate_final) == 0) {

    cat(
        "No duplicate ChEMBL assignments in final set.\n"
    )

} else {

    print(

        duplicate_final %>%

            select(
                final_rank,
                pert_id,
                pert_iname,
                chembl_id,
                final_evidence_score
            ),

        row.names = FALSE
    )

}


# ------------------------------------------------------------
# 19. SAVE FINAL RANKING
# ------------------------------------------------------------

write_csv(
    final_ranked,
    file.path(
        out_dir,
        "HER2_Luminal_FINAL_CANDIDATE_RANKING.csv"
    )
)


# ------------------------------------------------------------
# 20. SAVE TOP CANDIDATES
# ------------------------------------------------------------

write_csv(

    final_ranked %>%
        head(30),

    file.path(
        out_dir,
        "HER2_Luminal_TOP30_FINAL_CANDIDATES.csv"
    )
)


# ------------------------------------------------------------
# 21. SAVE PRIORITY A
# ------------------------------------------------------------

write_csv(

    priority_A,

    file.path(
        out_dir,
        "HER2_Luminal_PRIORITY_A_CANDIDATES.csv"
    )
)


# ------------------------------------------------------------
# 22. SAVE IDENTITY REVIEW SET
# ------------------------------------------------------------

write_csv(

    identity_review,

    file.path(
        out_dir,
        "HER2_Luminal_IDENTITY_REVIEW_CANDIDATES.csv"
    )
)


# ------------------------------------------------------------
# 23. SAVE NO-MECHANISM SET
# ------------------------------------------------------------

write_csv(

    strong_no_mechanism,

    file.path(
        out_dir,
        "HER2_Luminal_STRONG_LINCS_NO_MECHANISM.csv"
    )
)


# ------------------------------------------------------------
# 24. FINAL SUMMARY
# ------------------------------------------------------------

final_summary <- data.frame(

    total_candidates =
        nrow(final_ranked),

    priority_A =
        sum(
            final_ranked$final_priority ==
                "PRIORITY_A"
        ),

    priority_B =
        sum(
            final_ranked$final_priority ==
                "PRIORITY_B"
        ),

    priority_C =
        sum(
            final_ranked$final_priority ==
                "PRIORITY_C"
        ),

    exploratory =
        sum(
            final_ranked$final_priority ==
                "EXPLORATORY"
        ),

    high_confidence_identity =
        sum(
            final_ranked$identity_confidence ==
                "HIGH"
        ),

    identity_review =
        sum(
            final_ranked$identity_confidence ==
                "REVIEW"
        ),

    mechanism_supported =
        sum(
            final_ranked$
                mechanism_status_corrected ==
                "annotated"
        ),

    no_mechanism_record =
        sum(
            final_ranked$
                mechanism_status_corrected ==
                "no_mechanism_record"
        )

)

write_csv(
    final_summary,
    file.path(
        out_dir,
        "SCRIPT33_final_ranking_summary.csv"
    )
)


# ------------------------------------------------------------
# 25. COMPLETE
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 33 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Total candidates:",
    nrow(final_ranked),
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
    final_summary$high_confidence_identity,
    "\n"
)

cat(
    "Identity review required:",
    final_summary$identity_review,
    "\n"
)

cat(
    "Mechanism-supported candidates:",
    final_summary$mechanism_supported,
    "\n"
)

cat(
    "No mechanism record:",
    final_summary$no_mechanism_record,
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat("KEY OUTPUTS:\n\n")

cat(
    "1. HER2_Luminal_FINAL_CANDIDATE_RANKING.csv\n"
)

cat(
    "2. HER2_Luminal_TOP30_FINAL_CANDIDATES.csv\n"
)

cat(
    "3. HER2_Luminal_PRIORITY_A_CANDIDATES.csv\n"
)

cat(
    "4. HER2_Luminal_IDENTITY_REVIEW_CANDIDATES.csv\n"
)

cat(
    "5. HER2_Luminal_STRONG_LINCS_NO_MECHANISM.csv\n"
)

cat(
    "6. SCRIPT33_final_ranking_summary.csv\n\n"
)

cat(
    "IMPORTANT:\n"
)

cat(
    "This ranking represents computational drug-repurposing\n"
)

cat(
    "prioritization based on transcriptional reversal and\n"
)

cat(
    "supporting annotation evidence.\n"
)

cat(
    "It does NOT establish therapeutic efficacy.\n\n"
)

cat(
    "PIPELINE COMPLETE.\n"
)