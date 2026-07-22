# ============================================================
# VIOLIN PLOT MODULE
# ============================================================

#' Prepare data for violin plot
#' 
#' @param itx Intensity matrix (can be DelayedArray, matrix, or data.frame). 
#'           Rows are features, columns are pixels.
#' @param cluster_df Data frame with columns 'pixel_id' and the cluster annotation column.
#' @param feature Character, name of the metabolite/feature to plot.
#' @param cluster_col Character, name of the column in cluster_df containing group labels.
#' @return Data frame with columns: pixel_id, cluster, intensity.
violin_plot_data <- function(itx, cluster_df, feature, cluster_col) {
    
    pixel_ids <- cluster_df$pixel_id
    
    # Extract intensity values for the given feature and matching pixels
    # Handle DelayedArray, matrix, or data.frame
    if (inherits(itx, "DelayedArray")) {
        vals <- as.numeric(itx[feature, pixel_ids, drop = TRUE])
    } else if (is.matrix(itx) || is.data.frame(itx)) {
        vals <- as.numeric(itx[feature, pixel_ids])
    } else {
        stop("itx must be a matrix, data.frame, or DelayedArray")
    }
    
    # Build data frame
    df <- data.frame(
        pixel_id = pixel_ids,
        cluster = cluster_df[[cluster_col]],
        intensity = vals,
        stringsAsFactors = FALSE
    )
    
    # Remove any NA intensities (optional, but good practice)
    df <- df[!is.na(df$intensity), ]
    
    df
}

#' Create violin plot with overlaid boxplot
#' 
#' @param df Data frame from violin_plot_data().
#' @param feature Character, name of the feature (used for title).
#' @param x_label Character, label for x-axis (default "Cluster").
#' @param y_label Character, label for y-axis (default "Intensity").
#' @return ggplot object.
plot_violin <- function(df, feature, x_label = "Cluster", y_label = "Intensity") {
    
    ggplot(df, aes(x = cluster, y = intensity, fill = cluster)) +
        geom_violin(trim = TRUE) +
        geom_boxplot(width = 0.1, outlier.size = 0.3) +
        theme_bw() +
        labs(
            title = feature,
            x = x_label,
            y = y_label
        ) +
        theme(legend.position = "none")
}