# ============================================================
# SCRIPT 31C — LINCS ROBUST COMPOUND PRIORITIZATION
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Purpose:
#   Prioritize compounds using reproducibility across all
#   five breast cancer cell lines.
#
# Input:
#   Script 31 connectivity results
#
# Output:
#   Compound-level robustness ranking
#   Tiered candidate lists
#   QC plots
#
# IMPORTANT:
#   This script does NOT access the large LINCS GCTX file.
# ============================================================

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
})

cat("==============================================\n")
cat("SCRIPT 31C — LINCS ROBUST COMPOUND PRIORITIZATION\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

input_file <- paste0(
    "results/drug_repurposing/",
    "LINCS_HER2_Luminal/",
    "HER2_Luminal_LINCS_cell_line_connectivity.csv"
)

out_dir <- paste0(
    "results/drug_repurposing/",
    "LINCS_HER2_Luminal/"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. LOAD CONNECTIVITY RESULTS
# ------------------------------------------------------------

if (!file.exists(input_file)) {

    stop(
        "Script 31 connectivity file not found:\n",
        input_file
    )

}

cat("Loading Script 31 connectivity results...\n\n")

dat <- read_csv(
    input_file,
    show_col_types = FALSE
)

cat(
    "Rows:",
    nrow(dat),
    "\n"
)

cat(
    "Columns:",
    ncol(dat),
    "\n\n"
)


# ------------------------------------------------------------
# 3. REQUIRED COLUMN CHECK
# ------------------------------------------------------------

required_cols <- c(
    "sig_id",
    "pert_id",
    "pert_iname",
    "cell_id",
    "connectivity_score"
)

missing_cols <- setdiff(
    required_cols,
    colnames(dat)
)

if (length(missing_cols) > 0) {

    stop(
        "Missing required columns:\n",
        paste(missing_cols, collapse = ", ")
    )

}


# ------------------------------------------------------------
# 4. BASIC QC
# ------------------------------------------------------------

cat("==============================================\n")
cat("INPUT QC\n")
cat("==============================================\n\n")

cat(
    "Signatures:",
    nrow(dat),
    "\n"
)

cat(
    "Unique compounds:",
    n_distinct(dat$pert_id),
    "\n"
)

cat(
    "Cell lines:",
    n_distinct(dat$cell_id),
    "\n"
)

cat(
    "Missing connectivity:",
    sum(is.na(dat$connectivity_score)),
    "\n\n"
)

cell_lines <- sort(
    unique(dat$cell_id)
)

cat(
    "Cell lines represented:\n",
    paste(cell_lines, collapse = ", "),
    "\n\n"
)


# ------------------------------------------------------------
# 5. COMPOUND × CELL-LINE SUMMARY
# ------------------------------------------------------------

cat("==============================================\n")
cat("COMPOUND-LEVEL CONNECTIVITY SUMMARY\n")
cat("==============================================\n\n")

compound_summary <- dat %>%

    group_by(
        pert_id,
        pert_iname
    ) %>%

    summarise(

        cell_lines_tested =
            n_distinct(cell_id),

        positive_cell_lines =
            n_distinct(
                cell_id[connectivity_score > 0]
            ),

        negative_cell_lines =
            n_distinct(
                cell_id[connectivity_score < 0]
            ),

        median_score =
            median(
                connectivity_score,
                na.rm = TRUE
            ),

        mean_score =
            mean(
                connectivity_score,
                na.rm = TRUE
            ),

        min_score =
            min(
                connectivity_score,
                na.rm = TRUE
            ),

        max_score =
            max(
                connectivity_score,
                na.rm = TRUE
            ),

        sd_score =
            sd(
                connectivity_score,
                na.rm = TRUE
            ),

        .groups = "drop"
    ) %>%

    mutate(

        positive_fraction =
            positive_cell_lines /
            cell_lines_tested,

        negative_fraction =
            negative_cell_lines /
            cell_lines_tested

    )


# ------------------------------------------------------------
# 6. CELL-LINE COUNT QC
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("CELL-LINE COUNT QC\n")
cat("==============================================\n\n")

cat(
    "Maximum cell lines tested:",
    max(compound_summary$cell_lines_tested),
    "\n"
)

cat(
    "Maximum positive cell lines:",
    max(compound_summary$positive_cell_lines),
    "\n"
)

cat(
    "Maximum negative cell lines:",
    max(compound_summary$negative_cell_lines),
    "\n\n"
)

if (
    max(compound_summary$positive_cell_lines) >
    max(compound_summary$cell_lines_tested)
) {
    stop(
        "QC FAILED: positive cell-line count exceeds tested cell-line count."
    )
}

if (
    max(compound_summary$negative_cell_lines) >
    max(compound_summary$cell_lines_tested)
) {
    stop(
        "QC FAILED: negative cell-line count exceeds tested cell-line count."
    )
}

cat("Cell-line count QC PASSED.\n")


# ------------------------------------------------------------
# 7. POSITIVE FRACTION
# ------------------------------------------------------------

compound_summary <- compound_summary %>%

    mutate(

        positive_fraction =
            positive_cell_lines /
            cell_lines_tested,

        negative_fraction =
            negative_cell_lines /
            cell_lines_tested

    )


# ------------------------------------------------------------
# 8. FOCUS ON ALL-FIVE-CELL-LINE COMPOUNDS
# ------------------------------------------------------------

all_five <- compound_summary %>%

    filter(
        cell_lines_tested == length(cell_lines)
    )

cat(
    "Compounds tested in all",
    length(cell_lines),
    "cell lines:",
    nrow(all_five),
    "\n\n"
)


if (nrow(all_five) == 0) {

    stop(
        "No compounds were tested in all five cell lines."
    )

}


# ------------------------------------------------------------
# 9. ROBUST REVERSAL CRITERIA
# ------------------------------------------------------------

# Primary reproducibility definition:
#
#   - tested in all five cell lines
#   - positive connectivity in >= 3/5
#   - positive median connectivity
#
# This deliberately avoids requiring every cell line
# to be positive because biological heterogeneity is expected.

robust <- all_five %>%

    filter(
        positive_cell_lines >= 3,
        median_score > 0
    )


cat("==============================================\n")
cat("ROBUST REVERSAL SET\n")
cat("==============================================\n\n")

cat(
    "All-five compounds:",
    nrow(all_five),
    "\n"
)

cat(
    "Robust candidates:",
    nrow(robust),
    "\n\n"
)


# ------------------------------------------------------------
# 10. REPRODUCIBILITY SCORE
# ------------------------------------------------------------

# Reproducibility score:
#
#   positive_fraction × median connectivity
#
# This rewards both:
#   1. strength of reversal
#   2. consistency across cell lines

robust <- robust %>%

    mutate(

        reproducibility_score =
            positive_fraction *
            median_score

    )


# ------------------------------------------------------------
# 11. ROBUSTNESS TIER
# ------------------------------------------------------------

robust <- robust %>%

    mutate(

        tier = case_when(

            positive_cell_lines == 5 &
                min_score > 0 ~
                "Tier 1 — Consistent 5/5",

            positive_cell_lines >= 4 &
                median_score > 0 ~
                "Tier 2 — Strong 4/5+",

            positive_cell_lines >= 3 &
                median_score > 0 ~
                "Tier 3 — Moderate 3/5+",

            TRUE ~
                "Other"

        )

    )


# ------------------------------------------------------------
# 12. FINAL RANKING
# ------------------------------------------------------------

robust_ranked <- robust %>%

    arrange(
        desc(reproducibility_score),
        desc(median_score),
        desc(positive_cell_lines),
        desc(min_score)
    ) %>%

    mutate(
        robust_rank = row_number()
    )


# ------------------------------------------------------------
# 13. PRINT TOP CANDIDATES
# ------------------------------------------------------------

cat("==============================================\n")
cat("TOP ROBUST CANDIDATES\n")
cat("==============================================\n\n")

print(

    robust_ranked %>%

        select(
            robust_rank,
            pert_id,
            pert_iname,
            cell_lines_tested,
            positive_cell_lines,
            negative_cell_lines,
            median_score,
            mean_score,
            min_score,
            reproducibility_score,
            tier
        ) %>%

        head(30),

    n = 30
)


# ------------------------------------------------------------
# 14. TIER DISTRIBUTION
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TIER DISTRIBUTION\n")
cat("==============================================\n\n")

print(
    table(
        robust_ranked$tier
    )
)


# ------------------------------------------------------------
# 15. MOST CONSISTENT REVERSALS
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("CONSISTENT 5/5 REVERSALS\n")
cat("==============================================\n\n")

consistent_5of5 <- robust_ranked %>%

    filter(
        positive_cell_lines == 5,
        min_score > 0
    ) %>%

    arrange(
        desc(median_score)
    )

cat(
    "Candidates with positive connectivity in all five:",
    nrow(consistent_5of5),
    "\n\n"
)

print(

    consistent_5of5 %>%

        select(
            robust_rank,
            pert_id,
            pert_iname,
            positive_cell_lines,
            median_score,
            mean_score,
            min_score,
            reproducibility_score,
            tier
        ) %>%

        head(30),

    n = 30
)


# ------------------------------------------------------------
# 16. STRONGEST MEDIAN CONNECTIVITY
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("STRONGEST MEDIAN CONNECTIVITY\n")
cat("==============================================\n\n")

strongest <- robust_ranked %>%

    arrange(
        desc(median_score)
    )

print(

    strongest %>%

        select(
            robust_rank,
            pert_id,
            pert_iname,
            positive_cell_lines,
            median_score,
            mean_score,
            min_score,
            reproducibility_score,
            tier
        ) %>%

        head(30),

    n = 30
)


# ------------------------------------------------------------
# 17. SAVE ALL-FIVE COMPOUNDS
# ------------------------------------------------------------

write_csv(

    all_five,

    file.path(
        out_dir,
        "LINCS_all_5_cellline_compounds.csv"
    )

)


# ------------------------------------------------------------
# 18. SAVE ROBUST CANDIDATES
# ------------------------------------------------------------

write_csv(

    robust_ranked,

    file.path(
        out_dir,
        "HER2_Luminal_LINCS_ROBUST_COMPOUNDS.csv"
    )

)


# ------------------------------------------------------------
# 19. SAVE TIER 1
# ------------------------------------------------------------

tier1 <- robust_ranked %>%

    filter(
        tier == "Tier 1 — Consistent 5/5"
    )

write_csv(

    tier1,

    file.path(
        out_dir,
        "HER2_Luminal_LINCS_TIER1.csv"
    )

)


# ------------------------------------------------------------
# 20. SAVE TIER 2
# ------------------------------------------------------------

tier2 <- robust_ranked %>%

    filter(
        tier == "Tier 2 — Strong 4/5+"
    )

write_csv(

    tier2,

    file.path(
        out_dir,
        "HER2_Luminal_LINCS_TIER2.csv"
    )

)


# ------------------------------------------------------------
# 21. SAVE TIER 3
# ------------------------------------------------------------

tier3 <- robust_ranked %>%

    filter(
        tier == "Tier 3 — Moderate 3/5+"
    )

write_csv(

    tier3,

    file.path(
        out_dir,
        "HER2_Luminal_LINCS_TIER3.csv"
    )

)


# ------------------------------------------------------------
# 22. QC PLOT — MEDIAN VS REPRODUCIBILITY
# ------------------------------------------------------------

p1 <- ggplot(

    robust_ranked,

    aes(
        x = median_score,
        y = reproducibility_score,
        size = positive_cell_lines
    )

) +

    geom_point(alpha = 0.7) +

    labs(

        title =
            "HER2+ Luminal LINCS Robust Compound Prioritization",

        x =
            "Median Connectivity Score",

        y =
            "Reproducibility Score",

        size =
            "Positive Cell Lines"

    ) +

    theme_minimal()


ggsave(

    filename = file.path(
        out_dir,
        "HER2_Luminal_LINCS_robust_prioritization_QC.png"
    ),

    plot = p1,

    width = 9,
    height = 7,
    dpi = 300

)


# ------------------------------------------------------------
# 23. QC PLOT — TOP 20
# ------------------------------------------------------------

top20 <- robust_ranked %>%

    head(20) %>%

    mutate(

        label =
            reorder(
                pert_iname,
                reproducibility_score
            )

    )

p2 <- ggplot(

    top20,

    aes(
        x = reproducibility_score,
        y = label
    )

) +

    geom_col() +

    labs(

        title =
            "Top 20 Robust LINCS Reversal Candidates",

        x =
            "Reproducibility Score",

        y =
            "Compound"

    ) +

    theme_minimal()


ggsave(

    filename = file.path(
        out_dir,
        "HER2_Luminal_LINCS_top20_candidates.png"
    ),

    plot = p2,

    width = 10,
    height = 8,
    dpi = 300

)


# ------------------------------------------------------------
# 24. SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(

    total_signatures =
        nrow(dat),

    total_compounds =
        n_distinct(dat$pert_id),

    total_cell_lines =
        n_distinct(dat$cell_id),

    all_five_cellline_compounds =
        nrow(all_five),

    robust_candidates =
        nrow(robust_ranked),

    tier1_candidates =
        nrow(tier1),

    tier2_candidates =
        nrow(tier2),

    tier3_candidates =
        nrow(tier3)

)


write_csv(

    summary_df,

    file.path(
        out_dir,
        "SCRIPT31C_robust_prioritization_summary.csv"
    )

)


# ------------------------------------------------------------
# 25. FINAL MESSAGE
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 31C COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Total signatures:",
    nrow(dat),
    "\n"
)

cat(
    "All-five-cell-line compounds:",
    nrow(all_five),
    "\n"
)

cat(
    "Robust candidates:",
    nrow(robust_ranked),
    "\n"
)

cat(
    "Tier 1:",
    nrow(tier1),
    "\n"
)

cat(
    "Tier 2:",
    nrow(tier2),
    "\n"
)

cat(
    "Tier 3:",
    nrow(tier3),
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "NEXT STEP:\n",
    "Script 32 — candidate mechanism / target annotation\n"
)