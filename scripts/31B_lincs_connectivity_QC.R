# ============================================================
# SCRIPT 31B — LINCS CONNECTIVITY RANKING QC
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

cat("==============================================\n")
cat("SCRIPT 31B — LINCS CONNECTIVITY RANKING QC\n")
cat("==============================================\n\n")

out_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

file <- file.path(
    out_dir,
    "HER2_Luminal_LINCS_cell_line_connectivity.csv"
)

if (!file.exists(file)) {
    stop("Script 31 output not found:\n", file)
}

dat <- read_csv(
    file,
    show_col_types = FALSE
)

cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

cat("Available columns:\n")
print(colnames(dat))

required_cols <- c(
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
# 1. AGGREGATE WITHIN COMPOUND × CELL LINE
# ------------------------------------------------------------

compound_cellline <- dat %>%
    group_by(
        pert_id,
        pert_iname,
        cell_id
    ) %>%
    summarise(
        signatures = n(),

        median_connectivity = median(
            connectivity_score,
            na.rm = TRUE
        ),

        mean_connectivity = mean(
            connectivity_score,
            na.rm = TRUE
        ),

        best_connectivity = max(
            connectivity_score,
            na.rm = TRUE
        ),

        worst_connectivity = min(
            connectivity_score,
            na.rm = TRUE
        ),

        .groups = "drop"
    )

cat("\n==============================================\n")
cat("CELL-LINE CONNECTIVITY DISTRIBUTION\n")
cat("==============================================\n\n")

distribution <- compound_cellline %>%
    group_by(
        pert_id,
        pert_iname
    ) %>%
    summarise(
        cell_lines_tested = n_distinct(cell_id),
        .groups = "drop"
    ) %>%
    count(
        cell_lines_tested,
        name = "compounds"
    ) %>%
    arrange(
        cell_lines_tested
    )

print(distribution)

# ------------------------------------------------------------
# 2. CROSS-CELL-LINE SUMMARY
# ------------------------------------------------------------

compound_summary <- compound_cellline %>%
    group_by(
        pert_id,
        pert_iname
    ) %>%
    summarise(
        cell_lines_tested = n_distinct(cell_id),

        median_score = median(
            median_connectivity,
            na.rm = TRUE
        ),

        mean_score = mean(
            median_connectivity,
            na.rm = TRUE
        ),

        positive_cell_lines = sum(
            median_connectivity > 0,
            na.rm = TRUE
        ),

        negative_cell_lines = sum(
            median_connectivity < 0,
            na.rm = TRUE
        ),

        .groups = "drop"
    )

# ------------------------------------------------------------
# 3. ≥2 CELL LINES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP COMPOUNDS — ≥2 CELL LINES\n")
cat("==============================================\n\n")

top2 <- compound_summary %>%
    filter(
        cell_lines_tested >= 2
    ) %>%
    arrange(
        desc(median_score)
    )

print(
    top2 %>%
        select(
            pert_id,
            pert_iname,
            cell_lines_tested,
            positive_cell_lines,
            median_score,
            mean_score
        ) %>%
        head(30),
    row.names = FALSE
)

# ------------------------------------------------------------
# 4. ≥3 CELL LINES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP COMPOUNDS — ≥3 CELL LINES\n")
cat("==============================================\n\n")

top3 <- compound_summary %>%
    filter(
        cell_lines_tested >= 3
    ) %>%
    arrange(
        desc(median_score)
    )

print(
    top3 %>%
        select(
            pert_id,
            pert_iname,
            cell_lines_tested,
            positive_cell_lines,
            median_score,
            mean_score
        ) %>%
        head(30),
    row.names = FALSE
)

# ------------------------------------------------------------
# 5. ALL 5 CELL LINES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP COMPOUNDS — ALL 5 CELL LINES\n")
cat("==============================================\n\n")

all5 <- compound_summary %>%
    filter(
        cell_lines_tested == 5
    ) %>%
    arrange(
        desc(median_score)
    )

cat(
    "Compounds tested in all 5:",
    nrow(all5),
    "\n\n"
)

print(
    all5 %>%
        select(
            pert_id,
            pert_iname,
            cell_lines_tested,
            positive_cell_lines,
            negative_cell_lines,
            median_score,
            mean_score
        ) %>%
        head(30),
    row.names = FALSE
)

# ------------------------------------------------------------
# 6. ROBUST REVERSAL CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("ROBUST REVERSAL CANDIDATES\n")
cat("==============================================\n\n")

robust <- compound_summary %>%
    filter(
        cell_lines_tested >= 3,
        positive_cell_lines >= 3
    ) %>%
    arrange(
        desc(median_score)
    )

cat(
    "Candidates meeting criteria:",
    nrow(robust),
    "\n\n"
)

print(
    robust %>%
        select(
            pert_id,
            pert_iname,
            cell_lines_tested,
            positive_cell_lines,
            negative_cell_lines,
            median_score,
            mean_score
        ) %>%
        head(30),
    row.names = FALSE
)

# ------------------------------------------------------------
# 7. SAVE
# ------------------------------------------------------------

write_csv(
    compound_cellline,
    file.path(
        out_dir,
        "HER2_Luminal_connectivity_QC_compound_cellline.csv"
    )
)

write_csv(
    compound_summary,
    file.path(
        out_dir,
        "HER2_Luminal_connectivity_QC_compound_summary.csv"
    )
)

write_csv(
    robust,
    file.path(
        out_dir,
        "HER2_Luminal_ROBUST_CONNECTIVITY_CANDIDATES.csv"
    )
)

cat("\n==============================================\n")
cat("SCRIPT 31B COMPLETE\n")
cat("==============================================\n")