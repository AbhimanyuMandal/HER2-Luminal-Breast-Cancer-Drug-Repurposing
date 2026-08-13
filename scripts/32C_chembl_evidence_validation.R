# ============================================================
# SCRIPT 32C — ChEMBL EVIDENCE VALIDATION
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Purpose:
#   Correct ChEMBL mechanism/evidence classification from
#   Script 32/32B without removing or re-ranking candidates.
#
# Key principles:
#   1. Preserve all 66 robust LINCS candidates.
#   2. Distinguish identity confidence from mechanism evidence.
#   3. Distinguish actual mechanism records from
#      no-mechanism-record entries.
#   4. Explicitly flag duplicate ChEMBL assignments.
# ============================================================

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

cat("==============================================\n")
cat("SCRIPT 32C — ChEMBL EVIDENCE VALIDATION\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

out_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

identity_file <- file.path(
    out_dir,
    "HER2_Luminal_ChEMBL_compound_identity.csv"
)

mechanism_file <- file.path(
    out_dir,
    "HER2_Luminal_ChEMBL_mechanisms.csv"
)

candidate_file <- file.path(
    out_dir,
    "HER2_Luminal_candidate_mechanism_annotation.csv"
)

robust_file <- file.path(
    out_dir,
    "HER2_Luminal_robust_candidates.csv"
)


# ------------------------------------------------------------
# 2. CHECK INPUTS
# ------------------------------------------------------------

required_files <- c(
    identity_file,
    mechanism_file,
    candidate_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing required Script 32 outputs:\n",
        paste(missing_files, collapse = "\n")
    )

}


# ------------------------------------------------------------
# 3. LOAD SCRIPT 32 OUTPUTS
# ------------------------------------------------------------

cat("Loading Script 32 outputs...\n\n")

identity <- read_csv(
    identity_file,
    show_col_types = FALSE
)

mechanisms <- read_csv(
    mechanism_file,
    show_col_types = FALSE
)

candidate_annotation <- read_csv(
    candidate_file,
    show_col_types = FALSE
)

cat(
    "Identity records:",
    nrow(identity),
    "\n"
)

cat(
    "Mechanism records:",
    nrow(mechanisms),
    "\n"
)

cat(
    "Candidate annotation records:",
    nrow(candidate_annotation),
    "\n\n"
)


# ------------------------------------------------------------
# 4. INPUT STRUCTURE QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("INPUT STRUCTURE QC\n")
cat("==============================================\n\n")

cat("Identity columns:\n")
print(colnames(identity))

cat("\nMechanism columns:\n")
print(colnames(mechanisms))

cat("\nCandidate annotation columns:\n")
print(colnames(candidate_annotation))

cat("\n")


# ------------------------------------------------------------
# 5. NORMALIZE KEY IDENTIFIERS
# ------------------------------------------------------------

identity <- identity %>%
    mutate(
        chembl_id = as.character(chembl_id),
        pert_id = as.character(pert_id),
        pert_iname = as.character(pert_iname),
        match_status = as.character(match_status)
    )

mechanisms <- mechanisms %>%
    mutate(
        chembl_id = as.character(chembl_id),
        target_chembl_id = as.character(target_chembl_id),
        mechanism_of_action = as.character(mechanism_of_action),
        action_type = as.character(action_type),
        mechanism_status = as.character(mechanism_status)
    )

candidate_annotation <- candidate_annotation %>%
    mutate(
        pert_id = as.character(pert_id),
        pert_iname = as.character(pert_iname),
        chembl_id = as.character(chembl_id)
    )


# ------------------------------------------------------------
# 6. DEFINE ACTUAL MECHANISM RECORDS
# ------------------------------------------------------------
#
# IMPORTANT:
#
# Script 32B showed that some rows in the mechanism table
# contain a ChEMBL ID but have:
#
#   mechanism_status = "no_mechanism_record"
#
# These MUST NOT be counted as mechanism evidence.
#
# An actual mechanism record requires:
#
#   mechanism_status == "annotated"
#
# and a non-missing target/mechanism.
# ------------------------------------------------------------

mechanism_evidence <- mechanisms %>%
    filter(
        mechanism_status == "annotated",
        !is.na(target_chembl_id),
        target_chembl_id != ""
    ) %>%
    distinct(
        chembl_id,
        target_chembl_id,
        mechanism_of_action,
        action_type,
        .keep_all = TRUE
    )

cat("Actual mechanism-target records:",
    nrow(mechanism_evidence),
    "\n\n"
)


# ------------------------------------------------------------
# 7. COMPOUND-LEVEL MECHANISM STATUS
# ------------------------------------------------------------

mechanism_by_compound <- mechanism_evidence %>%
    group_by(chembl_id) %>%
    summarise(
        mechanism_record_count = n(),
        unique_targets = n_distinct(target_chembl_id),
        target_chembl_ids = paste(
            unique(target_chembl_id),
            collapse = ";"
        ),
        mechanisms = paste(
            unique(
                mechanism_of_action[
                    !is.na(mechanism_of_action)
                ]
            ),
            collapse = " | "
        ),
        action_types = paste(
            unique(
                action_type[
                    !is.na(action_type)
                ]
            ),
            collapse = ";"
        ),
        mechanism_status_corrected = "annotated",
        .groups = "drop"
    )


# ------------------------------------------------------------
# 8. IDENTITY + CORRECTED MECHANISM EVIDENCE
# ------------------------------------------------------------

validated_identity <- identity %>%
    left_join(
        mechanism_by_compound,
        by = "chembl_id"
    ) %>%
    mutate(

        mechanism_status_corrected = if_else(
            is.na(mechanism_status_corrected),
            "no_mechanism_record",
            mechanism_status_corrected
        ),

        mechanism_record_count = if_else(
            is.na(mechanism_record_count),
            0L,
            mechanism_record_count
        ),

        unique_targets = if_else(
            is.na(unique_targets),
            0L,
            unique_targets
        )
    )


# ------------------------------------------------------------
# 9. IDENTITY CONFIDENCE
# ------------------------------------------------------------

validated_identity <- validated_identity %>%
    mutate(

        identity_confidence = case_when(

            match_status == "exact_name_match" ~
                "HIGH",

            match_status == "search_match_requires_review" ~
                "REVIEW",

            TRUE ~
                "UNKNOWN"
        )
    )


# ------------------------------------------------------------
# 10. DUPLICATE ChEMBL ASSIGNMENT QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("DUPLICATE ChEMBL ASSIGNMENT QC\n")
cat("==============================================\n\n")

duplicate_ids <- validated_identity %>%
    filter(
        !is.na(chembl_id),
        chembl_id != ""
    ) %>%
    count(
        chembl_id,
        name = "n_LINCS_compounds"
    ) %>%
    filter(
        n_LINCS_compounds > 1
    )

cat(
    "Duplicate ChEMBL assignments:",
    nrow(duplicate_ids),
    "\n\n"
)

if (nrow(duplicate_ids) > 0) {

    duplicate_detail <- validated_identity %>%
        semi_join(
            duplicate_ids,
            by = "chembl_id"
        ) %>%
        group_by(chembl_id) %>%
        summarise(
            n_LINCS_compounds = n(),
            LINCS_compounds = paste(
                unique(pert_iname),
                collapse = "; "
            ),
            LINCS_ids = paste(
                unique(pert_id),
                collapse = "; "
            ),
            .groups = "drop"
        )

    print(duplicate_detail)

} else {

    duplicate_detail <- data.frame()

    cat("No duplicate ChEMBL assignments detected.\n")

}

cat("\n")


# ------------------------------------------------------------
# 11. COMPOUND-LEVEL EVIDENCE TABLE
# ------------------------------------------------------------

compound_level <- validated_identity %>%
    select(
        pert_id,
        pert_iname,
        chembl_id,
        pref_name,
        max_phase,
        match_status,
        identity_confidence,
        mechanism_status_corrected,
        mechanism_record_count,
        unique_targets,
        target_chembl_ids,
        mechanisms,
        action_types
    ) %>%
    distinct()


# ------------------------------------------------------------
# 12. CORRECTED MECHANISM QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("CORRECTED MECHANISM STATUS\n")
cat("==============================================\n\n")

mechanism_status_qc <- compound_level %>%
    count(
        mechanism_status_corrected,
        name = "compounds"
    ) %>%
    arrange(
        desc(compounds)
    )

print(mechanism_status_qc)

cat("\n")


# ------------------------------------------------------------
# 13. IDENTITY CONFIDENCE QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("IDENTITY CONFIDENCE\n")
cat("==============================================\n\n")

identity_qc <- compound_level %>%
    count(
        identity_confidence,
        name = "compounds"
    ) %>%
    arrange(
        desc(compounds)
    )

print(identity_qc)

cat("\n")


# ------------------------------------------------------------
# 14. ACTUAL MECHANISM TARGET QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("ACTUAL MECHANISM TARGETS\n")
cat("==============================================\n\n")

target_qc <- compound_level %>%
    filter(
        mechanism_status_corrected == "annotated"
    ) %>%
    select(
        pert_id,
        pert_iname,
        chembl_id,
        identity_confidence,
        mechanism_record_count,
        unique_targets,
        target_chembl_ids,
        mechanisms,
        action_types
    ) %>%
    arrange(
        desc(mechanism_record_count)
    )

cat(
    "Candidates with actual mechanism evidence:",
    nrow(target_qc),
    "\n\n"
)

print(
    head(target_qc, 30),
    row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. IMPORTANT DUPLICATE WARNING
# ------------------------------------------------------------

if (nrow(duplicate_ids) > 0) {

    cat("==============================================\n")
    cat("IMPORTANT — DUPLICATE ID WARNING\n")
    cat("==============================================\n\n")

    cat(
        "One or more ChEMBL IDs are assigned to multiple LINCS compounds.\n"
    )

    cat(
        "These compounds will NOT be merged or removed.\n"
    )

    cat(
        "The ambiguity will be carried into Script 33.\n\n"
    )

}


# ------------------------------------------------------------
# 16. MERGE WITH CANDIDATE ANNOTATION
# ------------------------------------------------------------

final_annotation <- candidate_annotation %>%
    select(
        everything(),
        -any_of(
            c(
                "identity_confidence",
                "mechanism_status_corrected",
                "mechanism_record_count",
                "unique_targets",
                "target_chembl_ids",
                "mechanisms",
                "action_types"
            )
        )
    ) %>%
    left_join(
        compound_level %>%
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
            ),
        by = c(
            "pert_id",
            "chembl_id"
        )
    )


# ------------------------------------------------------------
# 17. FINAL CANDIDATE COUNT QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("FINAL CANDIDATE COUNT QC\n")
cat("==============================================\n\n")

cat(
    "Candidate annotation rows:",
    nrow(final_annotation),
    "\n"
)

cat(
    "Unique LINCS compounds:",
    n_distinct(final_annotation$pert_id),
    "\n"
)

cat(
    "Unique ChEMBL IDs:",
    n_distinct(
        final_annotation$chembl_id[
            !is.na(final_annotation$chembl_id)
        ]
    ),
    "\n"
)

cat("\n")


# ------------------------------------------------------------
# 18. SAVE OUTPUTS
# ------------------------------------------------------------

write_csv(
    compound_level,
    file.path(
        out_dir,
        "SCRIPT32C_compound_evidence_QC.csv"
    )
)

write_csv(
    final_annotation,
    file.path(
        out_dir,
        "SCRIPT32C_candidate_evidence_annotation.csv"
    )
)

write_csv(
    mechanism_evidence,
    file.path(
        out_dir,
        "SCRIPT32C_actual_mechanism_records.csv"
    )
)

write_csv(
    duplicate_detail,
    file.path(
        out_dir,
        "SCRIPT32C_duplicate_ChEMBL_assignments.csv"
    )
)

write_csv(
    mechanism_status_qc,
    file.path(
        out_dir,
        "SCRIPT32C_mechanism_status_summary.csv"
    )
)

write_csv(
    identity_qc,
    file.path(
        out_dir,
        "SCRIPT32C_identity_confidence_summary.csv"
    )
)


# ------------------------------------------------------------
# 19. SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(

    candidates = n_distinct(
        final_annotation$pert_id
    ),

    unique_chembl_compounds = n_distinct(
        final_annotation$chembl_id[
            !is.na(final_annotation$chembl_id)
        ]
    ),

    high_confidence_identity = sum(
        compound_level$identity_confidence == "HIGH",
        na.rm = TRUE
    ),

    identity_review_required = sum(
        compound_level$identity_confidence == "REVIEW",
        na.rm = TRUE
    ),

    compounds_with_actual_mechanism = sum(
        compound_level$mechanism_status_corrected ==
            "annotated",
        na.rm = TRUE
    ),

    compounds_without_mechanism_record = sum(
        compound_level$mechanism_status_corrected ==
            "no_mechanism_record",
        na.rm = TRUE
    ),

    duplicate_chembl_assignments = nrow(
        duplicate_ids
    ),

    actual_mechanism_target_records =
        nrow(mechanism_evidence)

)

write_csv(
    summary_df,
    file.path(
        out_dir,
        "SCRIPT32C_validation_summary.csv"
    )
)


# ------------------------------------------------------------
# 20. COMPLETION
# ------------------------------------------------------------

cat("==============================================\n")
cat("SCRIPT 32C COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Candidates:",
    summary_df$candidates,
    "\n"
)

cat(
    "Unique ChEMBL compounds:",
    summary_df$unique_chembl_compounds,
    "\n"
)

cat(
    "High-confidence identities:",
    summary_df$high_confidence_identity,
    "\n"
)

cat(
    "Identity review required:",
    summary_df$identity_review_required,
    "\n"
)

cat(
    "Actual mechanism evidence:",
    summary_df$compounds_with_actual_mechanism,
    "\n"
)

cat(
    "No mechanism record:",
    summary_df$compounds_without_mechanism_record,
    "\n"
)

cat(
    "Duplicate ChEMBL assignments:",
    summary_df$duplicate_chembl_assignments,
    "\n"
)

cat(
    "Actual mechanism-target records:",
    summary_df$actual_mechanism_target_records,
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "NEXT STEP:\n",
    "Script 33 — final candidate validation and ranking\n"
)