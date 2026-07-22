# ============================================================
# VOLCANO PLOT MODULE
# ============================================================

#' Prepare data for volcano plot
#' 
#' @param res Data frame with columns: log_fc, adj_p_value, test_statistic, and a feature column.
#' @param feature_col Name of column containing feature identifiers (default "metabolite").
#' @param logfc_cutoff Numeric threshold for absolute log fold change (optional, used for grouping).
#' @param fdr_cutoff Numeric threshold for adjusted p-value (optional, used for grouping).
#' @return Data frame with additional columns: sig_group, y_value (user-selectable later).
prepare_volcano_data <- function(res, feature_col = "metabolite",
                                 logfc_cutoff = 0.1, fdr_cutoff = 0.05) {
    
    # Ensure required columns exist
    required_cols <- c("log_fc", "adj_p_value", "test_statistic")
    missing <- setdiff(required_cols, colnames(res))
    if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
    }
    
    res$feature <- res[[feature_col]]
    
    # Define significance groups based on cutoffs
    res$sig_group <- dplyr::case_when(
        abs(res$log_fc) >= logfc_cutoff & res$adj_p_value < fdr_cutoff & res$log_fc > 0 ~ "Up",
        abs(res$log_fc) >= logfc_cutoff & res$adj_p_value < fdr_cutoff & res$log_fc < 0 ~ "Down",
        TRUE ~ "Not significant"
    )
    
    # Pre-compute alternative y-values (but can be overridden in plot function)
    res$t_stat <- abs(res$test_statistic)
    res$neg_log10_p <- -log10(res$p_value)  # p_value must exist; if missing, use adj_p_value
    if (!"p_value" %in% colnames(res) && "adj_p_value" %in% colnames(res)) {
        res$neg_log10_p <- -log10(res$adj_p_value)
    }
    
    res
}

#' Create volcano plot with flexible y-axis and labeling
#'
#' @param res Data frame prepared by prepare_volcano_data().
#' @param top_n Number of top features to label (by absolute test statistic).
#' @param extra_features Character vector of additional feature names to label.
#' @param logfc_cutoff Numeric threshold for log fold change (vertical lines).
#' @param fdr_cutoff Numeric threshold for FDR (horizontal line if y = neglogp).
#' @param y_axis Character: "t" for |t-statistic|, "neglogp" for -log10(p-value).
#' @param point_alpha Numeric point transparency.
#' @param point_size Numeric point size.
#' @param label_size Numeric label text size.
#' 
#' @return ggplot object.
plot_volcano <- function(res, top_n = 20, extra_features = NULL,
                         logfc_cutoff = 0.1, fdr_cutoff = 0.05,
                         y_axis = "t", point_alpha = 0.6, point_size = 1.5,
                         label_size = 3) {
    
    # Determine y-axis column
    if (y_axis == "t") {
        y_col <- "t_stat"
        y_label <- "|t-statistic|"
    } else if (y_axis == "neglogp") {
        y_col <- "neg_log10_p"
        y_label <- expression(-log[10]("p-value"))
    } else {
        stop("y_axis must be 't' or 'neglogp'")
    }
    
    # Recompute significance groups based on current cutoffs (if changed)
    res <- res %>%
        dplyr::mutate(
            sig_group = dplyr::case_when(
                abs(log_fc) >= logfc_cutoff & adj_p_value < fdr_cutoff & log_fc > 0 ~ "Up",
                abs(log_fc) >= logfc_cutoff & adj_p_value < fdr_cutoff & log_fc < 0 ~ "Down",
                TRUE ~ "Not significant"
            )
        )
    
    # Select top features by the chosen y-value
    top_features <- res %>%
        dplyr::filter(sig_group != "Not significant") %>%
        dplyr::arrange(desc(.data[[y_col]])) %>%
        head(top_n)
    
    if (!is.null(extra_features)) {
        extra <- res %>% dplyr::filter(feature %in% extra_features)
        top_features <- dplyr::bind_rows(top_features, extra) %>% distinct()
    }
    
    # Plot
    p <- ggplot(res, aes(x = log_fc, y = .data[[y_col]], color = sig_group)) +
        geom_point(alpha = point_alpha, size = point_size) +
        scale_color_manual(values = c(
            "Up" = "red",
            "Down" = "blue",
            "Not significant" = "grey70"
        )) +
        geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
                   linetype = "dashed", color = "grey40") +
        theme_bw() +
        labs(
            x = expression(log[2] ~ "Fold Change"),
            y = y_label,
            title = "Volcano plot"
        )
    
    # Add horizontal line only if using -log10(p-value)
    if (y_axis == "neglogp") {
        p <- p + geom_hline(yintercept = -log10(fdr_cutoff),
                            linetype = "dashed", color = "grey40")
    }
    
    # Add labels
    if (nrow(top_features) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = top_features,
            aes(label = feature),
            size = label_size,
            max.overlaps = Inf
        )
    }
    
    p
}