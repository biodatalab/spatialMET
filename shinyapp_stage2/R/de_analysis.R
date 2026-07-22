# ============================================================
# DIFFERENTIAL ANALYSIS MODULE
# ============================================================

#' Run differential abundance tests between two groups
#'
#' @param itx Intensity matrix (features x pixels). Can be matrix or DelayedArray.
#' @param cluster_df Data frame with pixel_id and annotation column.
#' @param cluster_col Column name in cluster_df containing group labels.
#' @param group1 Name of first group.
#' @param group2 Name of second group.
#' @param method Character: "wilcoxon", "ttest", "hellinger", "limma_spatial"
#' @param pval_adj Adjustment method.
#' @param progress_callback Optional function to report progress (for Shiny).
#' @return Data frame with columns: metabolite, cluster, avg_log_itx1, avg_log_itx2,
#'         log_fc, test_statistic, p_value, adj_p_value.
run_de_analysis <- function(itx, cluster_df, cluster_col,
                            group1, group2, method = "wilcoxon", pval_adj = "fdr",
                            progress_callback = NULL) {
    
    report_progress <- function(step, detail = NULL) {
        if (!is.null(progress_callback)) progress_callback(step, detail)
    }
    
    report_progress(0.1, "Aligning pixel IDs")
    
    # ---- Handling "All_Others" for limma ----
    if (method == "limma_spatial" && group2 == "All_Others") {
        # Create a new column that pools all non-group1 pixels into "All_Others"
        cluster_df <- cluster_df %>%
            dplyr::mutate(
                .pooled_group = ifelse(.[[cluster_col]] == group1, group1, "All_Others")
            )
        # Override cluster_col to use the pooled column
        cluster_col <- ".pooled_group"
        # Ensure group1 and group2 are valid levels
        group1 <- group1
        group2 <- "All_Others"
        # For other methods, we could also support "All_Others" but they already have run_de_vs_all
    }
    
    # Extract pixel IDs and annotations
    df <- data.frame(
        pixel_id = as.character(cluster_df$pixel_id),
        annots = as.character(cluster_df[[cluster_col]]),
        stringsAsFactors = FALSE
    )
    df <- df[df$annots %in% c(group1, group2), ]
    if (nrow(df) == 0) {
        stop("No pixels found for the selected groups.")
    }
    df$annots <- factor(df$annots, levels = c(group1, group2))
    
    # Match pixel IDs
    itx_pixels <- colnames(itx)
    if (is.null(itx_pixels)) {
        stop("Intensity matrix must have column names (pixel IDs).")
    }
    itx_pixels <- as.character(itx_pixels)
    common_pixels <- intersect(itx_pixels, df$pixel_id)
    if (length(common_pixels) == 0) {
        stop("No common pixel IDs between intensity matrix and cluster data.")
    }
    
    df <- df[df$pixel_id %in% common_pixels, ]
    itx_sub <- itx[, common_pixels, drop = FALSE]
    itx_sub <- itx_sub[, df$pixel_id, drop = FALSE]
    
    report_progress(0.3, "Converting data")
    
    if (inherits(itx_sub, "DelayedArray")) {
        itx_sub <- as.matrix(itx_sub)
    }
    
    if (method == "limma_spatial") {
        report_progress(0.4, "Running spatial limma")
        coords <- cluster_df[match(df$pixel_id, cluster_df$pixel_id), c("x_coord", "y_coord")]
        res <- run_limma_spatial(itx_sub, df$annots, coords, group1, group2)
    } else {
        report_progress(0.4, "Running simple tests")
        group1_idx <- which(df$annots == group1)
        group2_idx <- which(df$annots == group2)
        res <- da_tests(
            x = itx_sub,
            group1 = group1_idx,
            group2 = group2_idx,
            group1_lbl = group1,
            group2_lbl = group2,
            method = method,
            pval_adj = pval_adj
        )
    }
    
    report_progress(1, "Done")
    return(res)
}