# ============================================================
# SCRIPT 35 — EXTERNAL CANDIDATE VALIDATION
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Purpose:
# Independently validate the highest-ranked LINCS candidates
# using ChEMBL and Europe PMC evidence.
#
# This is NOT a clinical efficacy assessment.
# ============================================================

suppressPackageStartupMessages({

    library(dplyr)
    library(readr)
    library(httr2)
    library(jsonlite)
    library(stringr)
    library(tibble)

})

cat("==============================================\n")
cat("SCRIPT 35 — EXTERNAL CANDIDATE VALIDATION\n")
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

out_dir <- file.path(
    base_dir,
    "external_validation"
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

if (!file.exists(ranking_file)) {

    stop(
        "Final ranking not found:\n",
        ranking_file
    )

}

if (!file.exists(mechanism_file)) {

    warning(
        "Mechanism-supported candidate file not found.\n",
        "Continuing using final ranking only."
    )

}

cat("Required ranking file found.\n\n")


# ------------------------------------------------------------
# 3. LOAD FINAL RANKING
# ------------------------------------------------------------

ranking <- read_csv(
    ranking_file,
    show_col_types = FALSE
)

cat(
    "Final candidates:",
    nrow(ranking),
    "\n"
)


# ------------------------------------------------------------
# 4. SELECT TOP 20
# ------------------------------------------------------------

top_candidates <- ranking %>%

    arrange(final_rank) %>%

    slice_head(n = 20)

cat(
    "Candidates selected for external validation:",
    nrow(top_candidates),
    "\n\n"
)

print(
    top_candidates %>%
        select(
            final_rank,
            pert_id,
            pert_iname,
            tier,
            positive_cell_lines,
            median_score,
            reproducibility_score,
            final_priority
        ),
    n = 20
)


# ------------------------------------------------------------
# 5. EUROPE PMC QUERY FUNCTION
# ------------------------------------------------------------

query_europepmc <- function(
    compound,
    timeout_seconds = 30
) {

    # Search specifically for the compound together with
    # breast cancer / HER2 terminology.

    query <- paste0(
        '"',
        compound,
        '" AND ',
        '(breast cancer OR breast carcinoma OR HER2 OR ERBB2)'
    )

    url <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    result <- tryCatch({

        req <- request(url) %>%

            req_url_query(
                query = query,
                format = "json",
                pageSize = 100
            ) %>%

            req_timeout(timeout_seconds)

        response <- req_perform(req)

        resp_body_json(response)

    }, error = function(e) {

        NULL

    })

    if (is.null(result)) {

        return(
            tibble(
                literature_status = "QUERY_FAILED",
                publication_count = NA_integer_
            )
        )

    }

    count <- result$hitCount

    if (is.null(count)) {

        count <- NA_integer_

    }

    tibble(

        literature_status =
            ifelse(
                is.na(count),
                "UNKNOWN",
                ifelse(
                    count > 0,
                    "EVIDENCE_FOUND",
                    "NO_MATCH_FOUND"
                )
            ),

        publication_count = count

    )

}


# ------------------------------------------------------------
# 6. GENERAL ONCOLOGY QUERY
# ------------------------------------------------------------

query_oncology <- function(
    compound,
    timeout_seconds = 30
) {

    query <- paste0(
        '"',
        compound,
        '" AND ',
        '(cancer OR tumor OR tumour OR oncology)'
    )

    url <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    result <- tryCatch({

        req <- request(url) %>%

            req_url_query(
                query = query,
                format = "json",
                pageSize = 100
            ) %>%

            req_timeout(timeout_seconds)

        response <- req_perform(req)

        resp_body_json(response)

    }, error = function(e) {

        NULL

    })

    if (is.null(result)) {

        return(
            tibble(
                oncology_status = "QUERY_FAILED",
                oncology_publications = NA_integer_
            )
        )

    }

    count <- result$hitCount

    if (is.null(count)) {

        count <- NA_integer_

    }

    tibble(

        oncology_status =
            ifelse(
                is.na(count),
                "UNKNOWN",
                ifelse(
                    count > 0,
                    "EVIDENCE_FOUND",
                    "NO_MATCH_FOUND"
                )
            ),

        oncology_publications = count

    )

}


# ------------------------------------------------------------
# 7. VALIDATE EACH CANDIDATE
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("EXTERNAL LITERATURE VALIDATION\n")
cat("==============================================\n\n")

validation_results <- vector(
    "list",
    nrow(top_candidates)
)


for (i in seq_len(nrow(top_candidates))) {

    compound <- top_candidates$pert_iname[i]

    cat(
        "[",
        i,
        "/",
        nrow(top_candidates),
        "] ",
        compound,
        "\n",
        sep = ""
    )

    breast_result <- query_europepmc(
        compound
    )

    Sys.sleep(0.2)

    oncology_result <- query_oncology(
        compound
    )

    validation_results[[i]] <-

        bind_cols(
            tibble(
                final_rank =
                    top_candidates$final_rank[i],

                pert_id =
                    top_candidates$pert_id[i],

                pert_iname =
                    compound,

                tier =
                    top_candidates$tier[i],

                final_priority =
                    top_candidates$final_priority[i],

                positive_cell_lines =
                    top_candidates$positive_cell_lines[i],

                median_score =
                    top_candidates$median_score[i],

                reproducibility_score =
                    top_candidates$reproducibility_score[i],

                chembl_id =
                    top_candidates$chembl_id[i],

                pref_name =
                    top_candidates$pref_name[i],

                identity_confidence =
                    top_candidates$identity_confidence[i],

                mechanism_status =
                    top_candidates$mechanism_status_corrected[i]

            ),

            breast_result,

            oncology_result
        )

}


validation <- bind_rows(
    validation_results
)


# ------------------------------------------------------------
# 8. EVIDENCE CLASSIFICATION
# ------------------------------------------------------------

validation <- validation %>%

    mutate(

        external_evidence_class = case_when(

            literature_status ==
                "EVIDENCE_FOUND" &
            oncology_status ==
                "EVIDENCE_FOUND" ~
                "BREAST_CANCER_OR_HER2_EVIDENCE",

            literature_status ==
                "EVIDENCE_FOUND" &
            oncology_status !=
                "EVIDENCE_FOUND" ~
                "BREAST_CANCER_HER2_SIGNAL",

            literature_status ==
                "NO_MATCH_FOUND" &
            oncology_status ==
                "EVIDENCE_FOUND" ~
                "GENERAL_ONCOLOGY_EVIDENCE",

            literature_status ==
                "NO_MATCH_FOUND" &
            oncology_status ==
                "NO_MATCH_FOUND" ~
                "LIMITED_EXTERNAL_EVIDENCE",

            TRUE ~
                "REVIEW_REQUIRED"

        )

    )


# ------------------------------------------------------------
# 9. COMPUTATIONAL + EXTERNAL EVIDENCE CLASS
# ------------------------------------------------------------

validation <- validation %>%

    mutate(

        overall_evidence_class = case_when(

            final_priority == "PRIORITY_A" &
            external_evidence_class ==
                "BREAST_CANCER_OR_HER2_EVIDENCE" ~
                "STRONG_CANDIDATE",

            final_priority %in%
                c("PRIORITY_A", "PRIORITY_B") &
            external_evidence_class ==
                "BREAST_CANCER_HER2_SIGNAL" ~
                "PROMISING_CANDIDATE",

            final_priority %in%
                c("PRIORITY_A", "PRIORITY_B") &
            external_evidence_class ==
                "GENERAL_ONCOLOGY_EVIDENCE" ~
                "MECHANISTICALLY_INTERESTING",

            final_priority %in%
                c("PRIORITY_A", "PRIORITY_B") &
            external_evidence_class ==
                "LIMITED_EXTERNAL_EVIDENCE" ~
                "COMPUTATIONAL_HYPOTHESIS",

            TRUE ~
                "REQUIRES_REVIEW"

        )

    )


# ------------------------------------------------------------
# 10. SAVE COMPLETE VALIDATION TABLE
# ------------------------------------------------------------

write_csv(
    validation,
    file.path(
        out_dir,
        "HER2_Luminal_TOP20_EXTERNAL_VALIDATION.csv"
    )
)


# ------------------------------------------------------------
# 11. SUMMARY
# ------------------------------------------------------------

summary_df <- validation %>%

    count(
        external_evidence_class,
        name = "candidates"
    ) %>%

    arrange(
        desc(candidates)
    )

write_csv(
    summary_df,
    file.path(
        out_dir,
        "SCRIPT35_external_evidence_summary.csv"
    )
)


# ------------------------------------------------------------
# 12. TOP CANDIDATES
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("TOP CANDIDATES — EXTERNAL EVIDENCE\n")
cat("==============================================\n\n")

print(

    validation %>%

        select(
            final_rank,
            pert_iname,
            final_priority,
            positive_cell_lines,
            median_score,
            reproducibility_score,
            chembl_id,
            identity_confidence,
            mechanism_status,
            literature_status,
            publication_count,
            oncology_status,
            oncology_publications,
            external_evidence_class,
            overall_evidence_class
        ) %>%

        arrange(final_rank),

    n = 20
)


# ------------------------------------------------------------
# 13. FINAL QC
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 35 QC\n")
cat("==============================================\n\n")

cat(
    "Candidates evaluated:",
    nrow(validation),
    "\n"
)

cat(
    "Literature queries completed:",
    sum(
        validation$literature_status !=
            "QUERY_FAILED",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Breast cancer/HER2 evidence:",
    sum(
        validation$external_evidence_class ==
            "BREAST_CANCER_OR_HER2_EVIDENCE",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "General oncology evidence:",
    sum(
        validation$external_evidence_class ==
            "GENERAL_ONCOLOGY_EVIDENCE",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Limited external evidence:",
    sum(
        validation$external_evidence_class ==
            "LIMITED_EXTERNAL_EVIDENCE",
        na.rm = TRUE
    ),
    "\n\n"
)


cat("==============================================\n")
cat("SCRIPT 35 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat(
    "NEXT STEP:\n",
    "Script 36 — final portfolio candidate summary\n"
)