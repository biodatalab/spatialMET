# ============================================================
# Metabolite mapping helper for pathway analysis
# ============================================================

#' Load a metabolite‑pathway database from a CSV file
#'
#' Required columns: metabolite_id, metabolite_name, mass, pathway
#'
#' @param file path to CSV file
#' @return data frame
load_metabolite_db <- function(file) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    required <- c("metabolite_id", "metabolite_name", "mass", "pathway")
    if (!all(required %in% colnames(df))) {
        stop("Database must contain columns: ", paste(required, collapse = ", "))
    }
    df$mass <- as.numeric(df$mass)
    df
}

#' Match query metabolites to database by mass tolerance and/or exact name
#'
#' @param query_vec character vector of query identifiers (e.g., mz values or names)
#' @param db_df data frame from load_metabolite_db()
#' @param mass_tol numeric mass tolerance (Da). If NULL, mass matching is skipped.
#' @param name_match logical, whether to match by exact name (case‑insensitive)
#' @return data frame with columns: query, matched_metabolite_id, matched_name, pathway, match_type
match_metabolites <- function(query_vec, db_df, mass_tol = 0.01, name_match = TRUE) {
    result <- data.frame(query = character(),
                         matched_metabolite_id = character(),
                         matched_name = character(),
                         pathway = character(),
                         match_type = character(),
                         stringsAsFactors = FALSE)
    
    for (q in query_vec) {
        matched <- FALSE
        # Name match (exact, case‑insensitive)
        if (name_match) {
            exact <- which(tolower(db_df$metabolite_name) == tolower(q))
            if (length(exact) > 0) {
                for (i in exact) {
                    result <- rbind(result, data.frame(query = q,
                                                       matched_metabolite_id = db_df$metabolite_id[i],
                                                       matched_name = db_df$metabolite_name[i],
                                                       pathway = db_df$pathway[i],
                                                       match_type = "exact_name",
                                                       stringsAsFactors = FALSE))
                }
                matched <- TRUE
                next
            }
        }
        # Mass matching (if query can be converted to numeric)
        if (!matched && !is.null(mass_tol) && !is.na(suppressWarnings(as.numeric(q)))) {
            mass_q <- as.numeric(q)
            within_tol <- which(abs(db_df$mass - mass_q) <= mass_tol)
            if (length(within_tol) > 0) {
                for (i in within_tol) {
                    result <- rbind(result, data.frame(query = q,
                                                       matched_metabolite_id = db_df$metabolite_id[i],
                                                       matched_name = db_df$metabolite_name[i],
                                                       pathway = db_df$pathway[i],
                                                       match_type = "mass_tolerance",
                                                       stringsAsFactors = FALSE))
                }
                matched <- TRUE
            }
        }
        if (!matched) {
            result <- rbind(result, data.frame(query = q,
                                               matched_metabolite_id = NA,
                                               matched_name = NA,
                                               pathway = NA,
                                               match_type = "no_match",
                                               stringsAsFactors = FALSE))
        }
    }
    result
}

#' Convert a matching result into a pathway list suitable for run_ora()
#'
#' @param match_df output from match_metabolites()
#' @return named list with elements `query` (character vector) and `pathway_list`
prepare_pathway_for_ora <- function(match_df) {
    # Keep only rows that have a pathway
    matched <- match_df[!is.na(match_df$pathway), ]
    if (nrow(matched) == 0) return(NULL)
    # Unique query metabolites (the original identifiers)
    query_metabs <- unique(matched$query)
    # Build pathway list: each pathway -> list of query metabolites
    pathway_list <- split(matched$query, matched$pathway)
    pathway_list <- lapply(pathway_list, unique)
    list(query = query_metabs, pathway_list = pathway_list)
}