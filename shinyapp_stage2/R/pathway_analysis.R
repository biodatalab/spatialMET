# ============================================================
# PATHWAY ANALYSIS MODULE
# ============================================================

#' Convert metabolite‑pathway annotation data frame to a named list
#'
#' @param ann_df data frame with columns `metabolite` and `pathway`
#' @return named list, each element is a character vector of metabolites in that pathway
prepare_pathway_list <- function(ann_df) {
    required <- c("metabolite", "pathway")
    if (!all(required %in% colnames(ann_df))) {
        stop("Annotation data frame must contain columns: ", paste(required, collapse = ", "))
    }
    # Remove empty strings and NA
    ann_df <- ann_df[!is.na(ann_df$metabolite) & !is.na(ann_df$pathway) & 
                         ann_df$metabolite != "" & ann_df$pathway != "", ]
    pathways <- split(ann_df$metabolite, ann_df$pathway)
    # Remove duplicate metabolites within each pathway
    pathways <- lapply(pathways, unique)
    return(pathways)
}

#' Perform over‑representation analysis (hypergeometric test)
#'
#' @param query character vector of metabolite names in the query set (e.g., top features)
#' @param background character vector of all metabolite names (universe)
#' @param pathway_list named list from prepare_pathway_list()
#' @param min_size minimum number of metabolites in a pathway to consider
#' @param max_size maximum number of metabolites in a pathway to consider
#' @param p_adj_method method for p‑value adjustment (passed to p.adjust)
#' @return data frame with columns: Pathway, Pathway_size, Overlap_count,
#'         Overlap_metabolites, p_value, adj_p_value
run_ora <- function(query, background, pathway_list,
                    min_size = 2, max_size = 200, p_adj_method = "BH") {
    
    if (length(query) == 0) return(NULL)
    query_unique <- unique(query)
    background_unique <- unique(background)
    total <- length(background_unique)
    
    results <- data.frame()
    for (pname in names(pathway_list)) {
        pathway_metabs <- unique(pathway_list[[pname]])
        pathway_size <- length(pathway_metabs)
        if (pathway_size < min_size || pathway_size > max_size) next
        
        overlap <- intersect(pathway_metabs, query_unique)
        overlap_count <- length(overlap)
        if (overlap_count == 0) next
        
        # Hypergeometric test (phyper)
        # Probability of observing overlap_count or more by chance
        p_val <- stats::phyper(overlap_count - 1, pathway_size, total - pathway_size,
                               length(query_unique), lower.tail = FALSE)
        
        results <- rbind(results, data.frame(
            Pathway = pname,
            Pathway_size = pathway_size,
            Overlap_count = overlap_count,
            Overlap_metabolites = paste(overlap, collapse = ";"),
            p_value = p_val,
            stringsAsFactors = FALSE
        ))
    }
    if (nrow(results) == 0) return(NULL)
    results$adj_p_value <- stats::p.adjust(results$p_value, method = p_adj_method)
    results <- results[order(results$adj_p_value, -results$Overlap_count), ]
    rownames(results) <- NULL
    return(results)
}

#' Create a horizontal bar plot of top enriched pathways
#'
#' @param res_df data frame from run_ora()
#' @param n_top maximum number of pathways to display
#' @param p_cutoff adjusted p‑value cutoff for including pathways
#' @return ggplot object, or NULL if no pathways pass the cutoff
plot_pathway_results <- function(res_df, n_top = 20, p_cutoff = 0.05) {
    if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
    df <- res_df[res_df$adj_p_value <= p_cutoff, ]
    if (nrow(df) == 0) return(NULL)
    df <- df[1:min(n_top, nrow(df)), ]
    df$neg_log_p <- -log10(df$adj_p_value)
    # Order for bar plot (ascending to have the largest at top)
    df <- df[order(df$neg_log_p), ]
    df$Pathway <- factor(df$Pathway, levels = df$Pathway)
    
    ggplot2::ggplot(df, ggplot2::aes(x = .data$neg_log_p, y = .data$Pathway, fill = .data$neg_log_p)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue") +
        ggplot2::theme_bw() +
        ggplot2::labs(x = expression(-log[10]("adjusted p-value")), y = "", title = "Top enriched pathways") +
        ggplot2::theme(legend.position = "none")
}

#' Format results for DT::datatable (optional)
#'
#' @param res_df data frame from run_ora()
#' @param p_cutoff adjusted p‑value cutoff for displaying
#' @return data frame ready for DT
format_pathway_table <- function(res_df, p_cutoff = 0.05) {
    if (is.null(res_df)) return(NULL)
    df <- res_df[res_df$adj_p_value <= p_cutoff, ]
    if (nrow(df) == 0) return(NULL)
    df <- df[, c("Pathway", "Pathway_size", "Overlap_count", "p_value", "adj_p_value")]
    colnames(df) <- c("Pathway", "Size", "Overlap", "p-value", "Adj. p-value")
    return(df)
}