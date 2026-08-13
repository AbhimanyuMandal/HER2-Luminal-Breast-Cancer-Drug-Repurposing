# ============================================================
# SCRIPT 34 — BIOLOGICAL CANDIDATE INTERPRETATION
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Purpose:
#   Integrate final LINCS ranking, ChEMBL mechanism evidence,
#   disease signature and pathway information.
#
# IMPORTANT:
#   This script DOES NOT re-rank candidates.
#   It does NOT claim therapeutic efficacy.
#   It generates biological interpretation tables.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
})

cat("==============================================\n")
cat("SCRIPT 34 — BIOLOGICAL CANDIDATE INTERPRETATION\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

base_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

ranking_file <- file.path(
    base_dir,
    "HER2_Luminal_FINAL_CANDIDATE_RANKING_QC.csv"
)

mechanism_file <- file.path(
    base_dir,
    "HER2_Luminal_MECHANISM_SUPPORTED_CANDIDATES.csv"
)

core_file <- paste0(
    "results/final_disease_signature/",
    "HER2_Luminal/",
    "HER2_Luminal_CORE_SIGNATURE.csv"
)

extended_file <- paste0(
    "results/final_disease_signature/",
    "HER2_Luminal/",
    "HER2_Luminal_EXTENDED_SIGNATURE.csv"
)

go_file <- paste0(
    "results/final_disease_signature/",
    "HER2_Luminal/",
    "HER2_Luminal_signature_GO_pathway_membership.csv"
)

out_dir <- file.path(
    base_dir,
    "biological_interpretation"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. CHECK INPUTS
# ------------------------------------------------------------

cat("Checking required files...\n")

required_files <- c(
    ranking_file,
    mechanism_file,
    core_file,
    extended_file,
    go_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    stop(
        "Missing required files:\n",
        paste(
            missing_files,
            collapse = "\n"
        )
    )

}

cat("All required files found.\n\n")


# ------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------

cat("Loading final ranking...\n")

ranking <- read_csv(
    ranking_file,
    show_col_types = FALSE
)

cat(
    "Final candidates:",
    nrow(ranking),
    "\n"
)


cat("Loading mechanism evidence...\n")

mechanism <- read_csv(
    mechanism_file,
    show_col_types = FALSE
)

cat(
    "Mechanism-supported records:",
    nrow(mechanism),
    "\n"
)


cat("Loading core signature...\n")

core <- read_csv(
    core_file,
    show_col_types = FALSE
)

cat(
    "Core genes:",
    nrow(core),
    "\n"
)


cat("Loading extended signature...\n")

extended <- read_csv(
    extended_file,
    show_col_types = FALSE
)

cat(
    "Extended genes:",
    nrow(extended),
    "\n"
)


cat("Loading GO pathway membership...\n")

go_map <- read_csv(
    go_file,
    show_col_types = FALSE
)

cat(
    "GO membership rows:",
    nrow(go_map),
    "\n\n"
)


# ------------------------------------------------------------
# 4. CORE DISEASE BIOLOGY
# ------------------------------------------------------------

cat("==============================================\n")
cat("CORE DISEASE BIOLOGY\n")
cat("==============================================\n\n")

cat("Core signature genes:\n")

print(
    core %>%
        select(
            any_of(
                c(
                    "gene",
                    "logFC",
                    "FDR",
                    "direction"
                )
            )
        ),
    row.names = FALSE
)


# ------------------------------------------------------------
# 5. CONTEXTUAL MARKERS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("CONTEXTUAL HER2 / LUMINAL MARKERS\n")
cat("==============================================\n\n")

contextual_genes <- extended %>%
    filter(
        gene %in% c(
            "ERBB2",
            "GRB7",
            "FOXA1",
            "AR",
            "ESR1",
            "PGR"
        )
    )

if (nrow(contextual_genes) > 0) {

    print(
        contextual_genes %>%
            select(
                any_of(
                    c(
                        "gene",
                        "logFC",
                        "FDR",
                        "direction"
                    )
                )
            ),
        row.names = FALSE
    )

} else {

    cat(
        "No contextual markers found in extended signature.\n"
    )

}


# ------------------------------------------------------------
# 6. PATHWAY SUMMARY
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("DISEASE PATHWAY SUMMARY\n")
cat("==============================================\n\n")

if (
    all(
        c(
            "gene",
            "pathway",
            "direction"
        ) %in%
        colnames(go_map)
    )
) {

    pathway_summary <- go_map %>%

        group_by(
            direction,
            pathway
        ) %>%

        summarise(
            genes = n_distinct(gene),
            gene_list = paste(
                unique(gene),
                collapse = ";"
            ),
            .groups = "drop"
        ) %>%

        arrange(
            direction,
            desc(genes)
        )

    print(
        pathway_summary %>%
            head(50),
        row.names = FALSE
    )

    write_csv(
        pathway_summary,
        file.path(
            out_dir,
            "HER2_Luminal_disease_pathway_summary.csv"
        )
    )

} else {

    warning(
        "GO pathway membership does not contain expected columns."
    )

}


# ------------------------------------------------------------
# 7. MECHANISM EVIDENCE SUMMARY
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("MECHANISM EVIDENCE SUMMARY\n")
cat("==============================================\n\n")

mechanism_summary <- mechanism %>%

    filter(
        !is.na(pert_id)
    ) %>%

    group_by(
        pert_id,
        pert_iname
    ) %>%

    summarise(

        chembl_id = first(
            na.omit(chembl_id),
            default = NA_character_
        ),

        pref_name = first(
            na.omit(pref_name),
            default = NA_character_
        ),

        mechanisms = paste(
            unique(
                na.omit(mechanism_of_action)
            ),
            collapse = " | "
        ),

        action_types = paste(
            unique(
                na.omit(action_type)
            ),
            collapse = " | "
        ),

        target_ids = paste(
            unique(
                na.omit(target_chembl_id)
            ),
            collapse = ";"
        ),

        mechanism_records = n(),

        .groups = "drop"
    )


cat(
    "Unique compounds with mechanism information:",
    nrow(mechanism_summary),
    "\n\n"
)


# ------------------------------------------------------------
# 8. JOIN MECHANISM TO FINAL RANKING
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("JOINING MECHANISM EVIDENCE TO FINAL RANKING\n")
cat("==============================================\n\n")

# Rename mechanism columns before joining so that they can
# never collide with columns already present in ranking.

mechanism_summary_join <- mechanism_summary %>%

    rename(
        mech_chembl_id = chembl_id,
        mech_pref_name = pref_name,
        mech_mechanisms = mechanisms,
        mech_action_types = action_types,
        mech_target_ids = target_ids,
        mech_mechanism_records = mechanism_records
    )


candidate_interpretation <- ranking %>%

    left_join(
        mechanism_summary_join,
        by = c(
            "pert_id",
            "pert_iname"
        )
    ) %>%

    mutate(

        mechanism_supported = case_when(

            !is.na(mech_mechanisms) &
                mech_mechanisms != "" ~ TRUE,

            TRUE ~ FALSE

        ),

        biological_evidence_class =
            case_when(

                final_priority == "PRIORITY_A" &
                    mechanism_supported ~
                    "HIGH_PRIORITY_MECHANISM_SUPPORTED",

                final_priority == "PRIORITY_A" &
                    !mechanism_supported ~
                    "HIGH_PRIORITY_NO_MECHANISM",

                final_priority == "PRIORITY_B" &
                    mechanism_supported ~
                    "STRONG_MECHANISM_SUPPORTED",

                final_priority == "PRIORITY_B" &
                    !mechanism_supported ~
                    "STRONG_NO_MECHANISM",

                final_priority == "PRIORITY_C" ~
                    "SECONDARY_CANDIDATE",

                TRUE ~
                    "EXPLORATORY"

            )
    ) %>%

    arrange(
        final_rank
    )


cat(
    "Candidates after mechanism join:",
    nrow(candidate_interpretation),
    "\n"
)

cat(
    "Candidates with mechanism evidence:",
    sum(
        candidate_interpretation$mechanism_supported,
        na.rm = TRUE
    ),
    "\n"
)

# Standardize mechanism column names for downstream analysis

candidate_interpretation <- candidate_interpretation %>%

    mutate(
        chembl_id = mech_chembl_id,
        pref_name = mech_pref_name,
        mechanisms = mech_mechanisms,
        action_types = mech_action_types,
        target_ids = mech_target_ids
    )

    
# ------------------------------------------------------------
# 9. TOP CANDIDATE INTERPRETATION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP CANDIDATE BIOLOGICAL INTERPRETATION\n")
cat("==============================================\n\n")

top_candidates <- candidate_interpretation %>%

    filter(
        final_rank <= 15
    ) %>%

    select(
        final_rank,
        pert_iname,
        pert_id,
        robust_rank,
        tier,
        positive_cell_lines,
        median_score,
        mean_score,
        reproducibility_score,
        identity_confidence,
        chembl_id,
        pref_name,
        mechanisms,
        action_types,
        biological_evidence_class
    )

print(
    top_candidates,
    row.names = FALSE
)


# ------------------------------------------------------------
# 10. PRIORITY A INTERPRETATION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("PRIORITY A BIOLOGICAL INTERPRETATION\n")
cat("==============================================\n\n")

priority_a <- candidate_interpretation %>%

    filter(
        final_priority == "PRIORITY_A"
    ) %>%

    select(
        final_rank,
        pert_iname,
        pert_id,
        positive_cell_lines,
        median_score,
        mean_score,
        reproducibility_score,
        identity_confidence,
        chembl_id,
        pref_name,
        mechanisms,
        action_types,
        target_ids,
        biological_evidence_class
    )

print(
    priority_a,
    row.names = FALSE
)


# ------------------------------------------------------------
# 11. MECHANISM-SUPPORTED TOP CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP MECHANISM-SUPPORTED CANDIDATES\n")
cat("==============================================\n\n")

mechanism_supported_top <- candidate_interpretation %>%

    filter(
        mechanism_supported
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
        pref_name,
        mechanisms,
        action_types,
        target_ids
    )

print(
    mechanism_supported_top %>%
        head(30),
    row.names = FALSE
)


# ------------------------------------------------------------
# 12. MECHANISM CATEGORY CLASSIFICATION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("MECHANISM CATEGORY SUMMARY\n")
cat("==============================================\n\n")

candidate_interpretation <- candidate_interpretation %>%

    mutate(

        mechanism_category = case_when(

            str_detect(
                str_to_lower(
                    mechanisms
                ),
                "kinase"
            ) ~ "KINASE",

            str_detect(
                str_to_lower(
                    mechanisms
                ),
                "parp"
            ) ~ "DNA_DAMAGE_RESPONSE",

            str_detect(
                str_to_lower(
                    mechanisms
                ),
                "tubulin|microtubule"
            ) ~ "CYTOSKELETON",

            str_detect(
                str_to_lower(
                    mechanisms
                ),
                "estrogen|androgen|hormone"
            ) ~ "HORMONE_SIGNALING",

            str_detect(
                str_to_lower(
                    mechanisms
                ),
                "vegf|vascular endothelial"
            ) ~ "ANGIOGENESIS",

            str_detect(
                str_to_lower(
                    mechanisms
                ),
                "mucin"
            ) ~ "CELL_SURFACE",

            TRUE ~ "OTHER / UNCLASSIFIED"

        )
    )


mechanism_category_summary <- candidate_interpretation %>%

    filter(
        mechanism_supported
    ) %>%

    count(
        mechanism_category,
        sort = TRUE
    )

print(
    mechanism_category_summary,
    row.names = FALSE
)


# ------------------------------------------------------------
# 13. FINAL INTERPRETATION TABLE
# ------------------------------------------------------------

interpretation_table <- candidate_interpretation %>%

    select(
        final_rank,
        pert_id,
        pert_iname,
        robust_rank,
        tier,
        positive_cell_lines,
        negative_cell_lines,
        median_score,
        mean_score,
        min_score,
        reproducibility_score,
        identity_confidence,
        chembl_id,
        pref_name,
        mechanisms,
        action_types,
        target_ids,
        mechanism_supported,
        mechanism_category,
        biological_evidence_class,
        final_evidence_score,
        final_priority
    )

write_csv(
    interpretation_table,
    file.path(
        out_dir,
        "HER2_Luminal_FINAL_BIOLOGICAL_INTERPRETATION.csv"
    )
)


# ------------------------------------------------------------
# 14. SAVE PRIORITY A
# ------------------------------------------------------------

write_csv(
    priority_a,
    file.path(
        out_dir,
        "HER2_Luminal_PRIORITY_A_BIOLOGICAL_INTERPRETATION.csv"
    )
)


# ------------------------------------------------------------
# 15. SAVE MECHANISM SUMMARY
# ------------------------------------------------------------

write_csv(
    mechanism_category_summary,
    file.path(
        out_dir,
        "HER2_Luminal_MECHANISM_CATEGORY_SUMMARY.csv"
    )
)


# ------------------------------------------------------------
# 16. SAVE TOP 15
# ------------------------------------------------------------

write_csv(
    top_candidates,
    file.path(
        out_dir,
        "HER2_Luminal_TOP15_BIOLOGICAL_INTERPRETATION.csv"
    )
)


# ------------------------------------------------------------
# 17. FINAL REPORT
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 34 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Final candidates:",
    nrow(candidate_interpretation),
    "\n"
)

cat(
    "Mechanism-supported:",
    sum(
        candidate_interpretation$mechanism_supported,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Priority A:",
    sum(
        candidate_interpretation$final_priority ==
            "PRIORITY_A",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Priority B:",
    sum(
        candidate_interpretation$final_priority ==
            "PRIORITY_B",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Priority C:",
    sum(
        candidate_interpretation$final_priority ==
            "PRIORITY_C",
        na.rm = TRUE
    ),
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat("KEY OUTPUTS:\n")
cat("1. HER2_Luminal_FINAL_BIOLOGICAL_INTERPRETATION.csv\n")
cat("2. HER2_Luminal_PRIORITY_A_BIOLOGICAL_INTERPRETATION.csv\n")
cat("3. HER2_Luminal_MECHANISM_CATEGORY_SUMMARY.csv\n")
cat("4. HER2_Luminal_TOP15_BIOLOGICAL_INTERPRETATION.csv\n")
cat("5. HER2_Luminal_disease_pathway_summary.csv\n\n")

cat("IMPORTANT:\n")
cat(
    "This script integrates computational reversal and\n",
    "mechanistic annotation. It does NOT establish clinical\n",
    "efficacy or therapeutic suitability.\n"
)

cat("\nNEXT STEP:\n")
cat(
    "Review the Priority A and Top 15 biological interpretation\n",
    "tables before drawing mechanistic conclusions.\n"
)