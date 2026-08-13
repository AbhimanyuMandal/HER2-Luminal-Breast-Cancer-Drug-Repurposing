# ============================================================
# SCRIPT 32 — CANDIDATE MECHANISM / TARGET ANNOTATION
# HER2+ Luminal Breast Cancer Drug Repurposing
#
# Purpose:
#   Annotate robust LINCS candidates with ChEMBL
#   compound and mechanism/target information.
#
# Important:
#   - Does NOT access the 22-GB GCTX file.
#   - Does NOT silently discard compounds without annotation.
#   - ChEMBL is used for annotation, not for final biological
#     prioritization.
# ============================================================

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(jsonlite)
})

cat("==============================================\n")
cat("SCRIPT 32 — CANDIDATE MECHANISM / TARGET ANNOTATION\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

base_dir <- "results/drug_repurposing/LINCS_HER2_Luminal"

candidate_file <- file.path(
    base_dir,
    "HER2_Luminal_LINCS_ROBUST_COMPOUNDS.csv"
)

pert_info_file <- file.path(
    "/Volumes/T7/LINCS/GSE92742",
    "GSE92742_Broad_LINCS_pert_info.txt.gz"
)

disease_signature_file <- paste0(
    "results/final_disease_signature/",
    "HER2_Luminal/",
    "HER2_Luminal_FINAL_RANKED_SIGNATURE.csv"
)

core_signature_file <- paste0(
    "results/final_disease_signature/",
    "HER2_Luminal/",
    "HER2_Luminal_CORE_SIGNATURE.csv"
)

out_dir <- base_dir

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. CHECK LOCAL FILES
# ------------------------------------------------------------

if (!file.exists(candidate_file)) {

    stop(
        "Robust candidate file not found:\n",
        candidate_file
    )

}

if (!file.exists(disease_signature_file)) {

    stop(
        "Final disease signature not found:\n",
        disease_signature_file
    )

}

if (!file.exists(core_signature_file)) {

    stop(
        "Core disease signature not found:\n",
        core_signature_file
    )

}


# ------------------------------------------------------------
# 3. LOAD CANDIDATES
# ------------------------------------------------------------

cat("Loading Script 31C robust candidates...\n")

candidates <- read_csv(
    candidate_file,
    show_col_types = FALSE
)

cat(
    "Robust candidates:",
    nrow(candidates),
    "\n\n"
)


required_candidate_cols <- c(
    "pert_id",
    "pert_iname",
    "cell_lines_tested",
    "positive_cell_lines",
    "median_score",
    "mean_score",
    "min_score",
    "reproducibility_score",
    "tier"
)

missing_cols <- setdiff(
    required_candidate_cols,
    colnames(candidates)
)

if (length(missing_cols) > 0) {

    stop(
        "Candidate file is missing:\n",
        paste(missing_cols, collapse = ", ")
    )

}


# ------------------------------------------------------------
# 4. LOAD DISEASE SIGNATURES
# ------------------------------------------------------------

cat("Loading disease signatures...\n")

disease_sig <- read_csv(
    disease_signature_file,
    show_col_types = FALSE
)

core_sig <- read_csv(
    core_signature_file,
    show_col_types = FALSE
)

cat(
    "Extended disease signature:",
    nrow(disease_sig),
    "genes\n"
)

cat(
    "Core signature:",
    nrow(core_sig),
    "genes\n\n"
)


# ------------------------------------------------------------
# 5. LOCAL PERTURBATION METADATA
# ------------------------------------------------------------

cat("Checking LINCS perturbation metadata...\n")

if (file.exists(pert_info_file)) {

    pert_info <- read_tsv(
        pert_info_file,
        show_col_types = FALSE
    )

    cat(
        "LINCS perturbation records:",
        nrow(pert_info),
        "\n\n"
    )

    candidates <- candidates %>%

        left_join(
            pert_info %>%
                select(
                    pert_id,
                    canonical_smiles,
                    inchi_key,
                    inchi_key_prefix,
                    pubchem_cid
                ) %>%
                distinct(pert_id, .keep_all = TRUE),
            by = "pert_id"
        )

} else {

    cat(
        "LINCS pert_info not found locally.\n"
    )

    cat(
        "Continuing without chemical identifiers.\n\n"
    )

    candidates$canonical_smiles <- NA_character_
    candidates$inchi_key <- NA_character_
    candidates$inchi_key_prefix <- NA_character_
    candidates$pubchem_cid <- NA_character_

}


# ------------------------------------------------------------
# 6. DISEASE SIGNATURE GENE SETS
# ------------------------------------------------------------

disease_genes <- unique(
    disease_sig$gene
)

core_genes <- unique(
    core_sig$gene
)

up_genes <- unique(
    disease_sig$gene[
        disease_sig$direction == "UP"
    ]
)

down_genes <- unique(
    disease_sig$gene[
        disease_sig$direction == "DOWN"
    ]
)

cat("Disease genes:", length(disease_genes), "\n")
cat("Core genes:", length(core_genes), "\n")
cat("Disease UP:", length(up_genes), "\n")
cat("Disease DOWN:", length(down_genes), "\n\n")


# ------------------------------------------------------------
# 7. CHOOSE API IDENTIFIER
# ------------------------------------------------------------

# Prefer chemical identifiers over compound names.
#
# ChEMBL can search molecules by name and also supports
# structure-based searches. We initially use the compound
# name for broad identity matching and retain chemical
# identifiers for future validation.

cat("==============================================\n")
cat("CHemBL COMPOUND IDENTIFICATION\n")
cat("==============================================\n\n")


# ------------------------------------------------------------
# 8. HELPER: SAFE JSON GET
# ------------------------------------------------------------

safe_get_json <- function(url) {

    result <- tryCatch(

        {

            txt <- readLines(
                url,
                warn = FALSE,
                encoding = "UTF-8"
            )

            txt <- paste(
                txt,
                collapse = "\n"
            )

            fromJSON(
                txt,
                simplifyVector = FALSE
            )

        },

        error = function(e) {

            NULL

        }

    )

    result

}


# ------------------------------------------------------------
# 9. HELPER: URL ENCODING
# ------------------------------------------------------------

url_encode <- function(x) {

    utils::URLencode(
        x,
        reserved = TRUE
    )

}


# ------------------------------------------------------------
# 10. CHemBL MOLECULE SEARCH
# ------------------------------------------------------------

search_chembl_molecule <- function(compound_name) {

    if (
        is.na(compound_name) ||
        compound_name == ""
    ) {

        return(
            data.frame(
                chembl_id = NA_character_,
                pref_name = NA_character_,
                max_phase = NA_real_,
                match_status = "missing_name",
                stringsAsFactors = FALSE
            )
        )

    }

    query_url <- paste0(
        "https://www.ebi.ac.uk/chembl/api/data/",
        "molecule/search.json?q=",
        url_encode(compound_name)
    )

    obj <- safe_get_json(
        query_url
    )

    if (is.null(obj)) {

        return(
            data.frame(
                chembl_id = NA_character_,
                pref_name = NA_character_,
                max_phase = NA_real_,
                match_status = "api_error",
                stringsAsFactors = FALSE
            )
        )

    }

    molecules <- obj$molecules

    if (
        is.null(molecules) ||
        length(molecules) == 0
    ) {

        return(
            data.frame(
                chembl_id = NA_character_,
                pref_name = NA_character_,
                max_phase = NA_real_,
                match_status = "no_match",
                stringsAsFactors = FALSE
            )
        )

    }

    # Extract first exact-ish candidate.
    #
    # We do NOT claim this is necessarily the correct
    # identity. It is flagged for downstream QC.

    ids <- vapply(
        molecules,
        function(x) {
            if (
                is.null(x$molecule_chembl_id)
            ) NA_character_
            else x$molecule_chembl_id
        },
        character(1)
    )

    names <- vapply(
        molecules,
        function(x) {
            if (
                is.null(x$pref_name)
            ) NA_character_
            else x$pref_name
        },
        character(1)
    )

    phases <- vapply(
        molecules,
        function(x) {

            if (
                is.null(x$max_phase)
            ) {

                NA_real_

            } else {

                suppressWarnings(
                    as.numeric(x$max_phase)
                )

            }

        },
        numeric(1)
    )

    exact_idx <- which(
        tolower(names) ==
        tolower(compound_name)
    )

    if (length(exact_idx) > 0) {

        idx <- exact_idx[1]
        status <- "exact_name_match"

    } else {

        idx <- 1
        status <- "search_match_requires_review"

    }

    data.frame(

        chembl_id =
            ids[idx],

        pref_name =
            names[idx],

        max_phase =
            phases[idx],

        match_status =
            status,

        stringsAsFactors =
            FALSE

    )

}


# ------------------------------------------------------------
# 11. RUN MOLECULE IDENTIFICATION
# ------------------------------------------------------------

cat(
    "Querying ChEMBL for",
    nrow(candidates),
    "candidate compounds...\n\n"
)

chembl_identity_list <- vector(
    "list",
    nrow(candidates)
)

for (i in seq_len(nrow(candidates))) {

    compound_name <- candidates$pert_iname[i]

    if (
        is.na(compound_name) ||
        compound_name == ""
    ) {

        compound_name <- candidates$pert_id[i]

    }

    cat(
        sprintf(
            "[%d/%d] %s\n",
            i,
            nrow(candidates),
            compound_name
        )
    )

    result <- search_chembl_molecule(
        compound_name
    )

    result$pert_id <-
        candidates$pert_id[i]

    result$pert_iname <-
        candidates$pert_iname[i]

    chembl_identity_list[[i]] <-
        result

    # Small delay to avoid hammering the API
    Sys.sleep(0.15)

}

chembl_identity <- bind_rows(
    chembl_identity_list
)


# ------------------------------------------------------------
# 12. SAVE IDENTITY MAPPING
# ------------------------------------------------------------

write_csv(

    chembl_identity,

    file.path(
        out_dir,
        "HER2_Luminal_ChEMBL_compound_identity.csv"
    )

)


# ------------------------------------------------------------
# 13. IDENTITY QC
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("ChEMBL IDENTITY QC\n")
cat("==============================================\n\n")

print(
    table(
        chembl_identity$match_status,
        useNA = "ifany"
    )
)

cat("\n")


# ------------------------------------------------------------
# 14. MERGE IDENTITIES
# ------------------------------------------------------------

candidates_annotated <- candidates %>%

    left_join(
        chembl_identity %>%
            select(
                pert_id,
                chembl_id,
                pref_name,
                max_phase,
                match_status
            ),
        by = "pert_id"
    )


# ------------------------------------------------------------
# 15. HELPER: CHemBL MECHANISM QUERY
# ------------------------------------------------------------

get_mechanisms <- function(chembl_id) {

    if (
        is.na(chembl_id) ||
        chembl_id == ""
    ) {

        return(
            data.frame(
                chembl_id = NA_character_,
                target_chembl_id = NA_character_,
                mechanism_of_action = NA_character_,
                action_type = NA_character_,
                mechanism_status = "no_chembl_id",
                stringsAsFactors = FALSE
            )
        )

    }

    query_url <- paste0(
        "https://www.ebi.ac.uk/chembl/api/data/",
        "mechanism.json?",
        "molecule_chembl_id=",
        url_encode(chembl_id)
    )

    obj <- safe_get_json(
        query_url
    )

    if (is.null(obj)) {

        return(
            data.frame(
                chembl_id = chembl_id,
                target_chembl_id = NA_character_,
                mechanism_of_action = NA_character_,
                action_type = NA_character_,
                mechanism_status = "api_error",
                stringsAsFactors = FALSE
            )
        )

    }

    mechanisms <- obj$mechanisms

    if (
        is.null(mechanisms) ||
        length(mechanisms) == 0
    ) {

        return(
            data.frame(
                chembl_id = chembl_id,
                target_chembl_id = NA_character_,
                mechanism_of_action = NA_character_,
                action_type = NA_character_,
                mechanism_status = "no_mechanism_record",
                stringsAsFactors = FALSE
            )
        )

    }

    rows <- lapply(

        mechanisms,

        function(m) {

            data.frame(

                chembl_id =
                    chembl_id,

                target_chembl_id =
                    if (
                        is.null(
                            m$target_chembl_id
                        )
                    )
                        NA_character_
                    else
                        m$target_chembl_id,

                mechanism_of_action =
                    if (
                        is.null(
                            m$mechanism_of_action
                        )
                    )
                        NA_character_
                    else
                        m$mechanism_of_action,

                action_type =
                    if (
                        is.null(
                            m$action_type
                        )
                    )
                        NA_character_
                    else
                        m$action_type,

                mechanism_status =
                    "annotated",

                stringsAsFactors =
                    FALSE

            )

        }

    )

    bind_rows(rows)

}


# ------------------------------------------------------------
# 16. QUERY MECHANISMS
# ------------------------------------------------------------

cat("==============================================\n")
cat("ChEMBL MECHANISM ANNOTATION\n")
cat("==============================================\n\n")

unique_chembl <- unique(
    na.omit(
        candidates_annotated$chembl_id
    )
)

cat(
    "Unique ChEMBL compounds:",
    length(unique_chembl),
    "\n\n"
)

mechanism_list <- vector(
    "list",
    length(unique_chembl)
)

for (i in seq_along(unique_chembl)) {

    id <- unique_chembl[i]

    cat(
        sprintf(
            "[%d/%d] %s\n",
            i,
            length(unique_chembl),
            id
        )
    )

    mechanism_list[[i]] <-
        get_mechanisms(id)

    Sys.sleep(0.15)

}

if (length(mechanism_list) > 0) {

    mechanism_table <-
        bind_rows(mechanism_list)

} else {

    mechanism_table <-
        data.frame()

}


# ------------------------------------------------------------
# 17. SAVE MECHANISMS
# ------------------------------------------------------------

write_csv(

    mechanism_table,

    file.path(
        out_dir,
        "HER2_Luminal_ChEMBL_mechanisms.csv"
    )

)


# ------------------------------------------------------------
# 18. TARGET OVERLAP
# ------------------------------------------------------------

# ChEMBL mechanism records provide target_chembl_id,
# but a target is not automatically equivalent to a gene
# symbol. We therefore retain target IDs and mechanism text
# here rather than inventing gene mappings.

target_summary <- mechanism_table %>%

    filter(
        !is.na(target_chembl_id)
    ) %>%

    distinct(
        chembl_id,
        target_chembl_id,
        mechanism_of_action,
        action_type
    )


# ------------------------------------------------------------
# 19. MERGE CANDIDATE + MECHANISM INFORMATION
# ------------------------------------------------------------

candidate_mechanisms <- candidates_annotated %>%

    left_join(
        target_summary,
        by = "chembl_id"
    )


# ------------------------------------------------------------
# 20. DISEASE SIGNATURE CONTEXT
# ------------------------------------------------------------

candidate_mechanisms <- candidate_mechanisms %>%

    mutate(

        disease_signature_genes =
            length(disease_genes),

        core_signature_genes =
            length(core_genes),

        disease_up_genes =
            length(up_genes),

        disease_down_genes =
            length(down_genes)

    )


# ------------------------------------------------------------
# 21. SAVE FULL ANNOTATION TABLE
# ------------------------------------------------------------

write_csv(

    candidate_mechanisms,

    file.path(
        out_dir,
        "HER2_Luminal_candidate_mechanism_annotation.csv"
    )

)


# ------------------------------------------------------------
# 22. ANNOTATION SUMMARY
# ------------------------------------------------------------

n_candidates <- nrow(
    candidates
)

n_chembl <- sum(
    !is.na(
        candidates_annotated$chembl_id
    )
)

n_exact <- sum(
    candidates_annotated$match_status ==
        "exact_name_match",
    na.rm = TRUE
)

n_mechanism <- n_distinct(
    mechanism_table$chembl_id[
        !is.na(
            mechanism_table$chembl_id
        )
    ]
)

n_target_records <- nrow(
    target_summary
)


summary_df <- data.frame(

    candidates =
        n_candidates,

    ChEMBL_compounds_identified =
        n_chembl,

    exact_name_matches =
        n_exact,

    compounds_with_mechanism_records =
        n_mechanism,

    mechanism_target_records =
        n_target_records,

    disease_signature_genes =
        length(disease_genes),

    core_signature_genes =
        length(core_genes)

)


write_csv(

    summary_df,

    file.path(
        out_dir,
        "SCRIPT32_annotation_summary.csv"
    )

)


# ------------------------------------------------------------
# 23. PRINT SUMMARY
# ------------------------------------------------------------

cat("\n==============================================\n")
cat("SCRIPT 32 COMPLETE\n")
cat("==============================================\n\n")

cat(
    "Robust candidates:",
    n_candidates,
    "\n"
)

cat(
    "ChEMBL compounds identified:",
    n_chembl,
    "\n"
)

cat(
    "Exact name matches:",
    n_exact,
    "\n"
)

cat(
    "Compounds with mechanism records:",
    n_mechanism,
    "\n"
)

cat(
    "Mechanism-target records:",
    n_target_records,
    "\n\n"
)

cat(
    "Output directory:\n",
    out_dir,
    "\n\n"
)

cat("Key outputs:\n\n")

cat(
    "1. HER2_Luminal_ChEMBL_compound_identity.csv\n"
)

cat(
    "2. HER2_Luminal_ChEMBL_mechanisms.csv\n"
)

cat(
    "3. HER2_Luminal_candidate_mechanism_annotation.csv\n"
)

cat(
    "4. SCRIPT32_annotation_summary.csv\n\n"
)

cat(
    "IMPORTANT:\n"
)

cat(
    "ChEMBL annotations are evidence for compound mechanism,\n"
)

cat(
    "not proof of HER2+ luminal therapeutic efficacy.\n\n"
)

cat(
    "NEXT STEP:\n"
)

cat(
    "Script 33 — final candidate validation and ranking\n"
)