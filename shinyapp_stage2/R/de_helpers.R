# R/de_helpers.R
# Parallel implementations of DE tests used in the app

library(BiocParallel)
library(matrixStats)

# Parallel backend (used by all functions below)
BPPARAM <- MulticoreParam(workers = max(1, availableCores() - 1))

#' Run differential analysis using Wilcoxon, t-test, or Hellinger (parallel)
#' @param itx intensity matrix (features x pixels)
#' @param cluster_df data frame with pixel_id and cluster column
#' @param cluster_col character: column name in cluster_df containing group labels
#' @param group1 character: first group label
#' @param group2 character: second group label
#' @param method "wilcoxon", "ttest", or "hellinger"
#' @param pval_adj p.adjust method (e.g., "fdr")
#' @param progress_callback function(step) called for progress
#' @return data frame with results
run_de_analysis_parallel <- function(itx, cluster_df, cluster_col, group1, group2, 
                                     method, pval_adj = "fdr", progress_callback = NULL) {
    # Prepare data
    df <- data.frame(pixel_id = cluster_df$pixel_id,
                     cluster = cluster_df[[cluster_col]],
                     stringsAsFactors = FALSE)
    df <- df[df$cluster %in% c(group1, group2), ]
    df$cluster <- factor(df$cluster, levels = c(group1, group2))
    group1_pixels <- df$pixel_id[df$cluster == group1]
    group2_pixels <- df$pixel_id[df$cluster == group2]
    
    itx_sub <- itx[, c(group1_pixels, group2_pixels), drop = FALSE]
    if (inherits(itx_sub, "DelayedArray")) {
        itx_sub <- as.matrix(itx_sub)
    }
    
    n1 <- length(group1_pixels)
    n2 <- length(group2_pixels)
    n_features <- nrow(itx_sub)
    
    # Pre-allocate
    log_fc <- numeric(n_features)
    avg1 <- numeric(n_features)
    avg2 <- numeric(n_features)
    p_vals <- numeric(n_features)
    stats <- numeric(n_features)
    
    # Parallel loop
    if (method == "wilcoxon") {
        results <- bplapply(seq_len(n_features), function(i) {
            x <- itx_sub[i, 1:n1]
            y <- itx_sub[i, (n1+1):(n1+n2)]
            test <- wilcox.test(x, y, exact = FALSE)
            logfc <- log2((median(x)+1e-6)/(median(y)+1e-6))
            list(logfc = logfc, avg1 = mean(x), avg2 = mean(y), 
                 stat = test$statistic, pval = test$p.value)
        }, BPPARAM = BPPARAM)
    } else if (method == "ttest") {
        results <- bplapply(seq_len(n_features), function(i) {
            x <- itx_sub[i, 1:n1]
            y <- itx_sub[i, (n1+1):(n1+n2)]
            test <- t.test(x, y)
            logfc <- log2((mean(x)+1e-6)/(mean(y)+1e-6))
            list(logfc = logfc, avg1 = mean(x), avg2 = mean(y), 
                 stat = test$statistic, pval = test$p.value)
        }, BPPARAM = BPPARAM)
    } else if (method == "hellinger") {
        results <- bplapply(seq_len(n_features), function(i) {
            x <- itx_sub[i, 1:n1]
            y <- itx_sub[i, (n1+1):(n1+n2)]
            logfc <- log2((mean(x)+1e-6)/(mean(y)+1e-6))
            # Approximate Hellinger distance between two normal distributions
            mu1 <- mean(x); sd1 <- sd(x)
            mu2 <- mean(y); sd2 <- sd(y)
            hell <- 1 - sqrt(2*sd1*sd2/(sd1^2+sd2^2)) * exp(-0.25*(mu1-mu2)^2/(sd1^2+sd2^2))
            # For Hellinger we don't have a proper p‑value; we'll set it to NA
            # and let the calling code handle it
            list(logfc = logfc, avg1 = mu1, avg2 = mu2, 
                 stat = hell, pval = NA_real_)
        }, BPPARAM = BPPARAM)
        warning("Hellinger test: p‑values not available, using NA.")
    } else {
        stop("Unknown test method")
    }
    
    # Extract results
    for (i in seq_len(n_features)) {
        res <- results[[i]]
        log_fc[i] <- res$logfc
        avg1[i] <- res$avg1
        avg2[i] <- res$avg2
        stats[i] <- res$stat
        p_vals[i] <- if (is.null(res$pval)) NA else res$pval
        if (!is.null(progress_callback)) progress_callback(i / n_features)
    }
    
    # Adjust p‑values if any are non‑NA
    if (any(!is.na(p_vals))) {
        adj_p <- p.adjust(p_vals, method = pval_adj)
    } else {
        adj_p <- rep(NA, n_features)
    }
    
    res_df <- data.frame(
        metabolite = rownames(itx_sub),
        cluster = paste(group1, "vs", group2),
        avg_log_itx1 = avg1,
        avg_log_itx2 = avg2,
        log_fc = log_fc,
        test_statistic = stats,
        p_value = p_vals,
        adj_p_value = adj_p,
        stringsAsFactors = FALSE
    )
    return(res_df)
}

#' Run DE of one group vs all others (parallel)
#' @param ... passed to run_de_analysis_parallel
#' @return data frame with results, cluster column modified
run_de_vs_all <- function(itx, cluster_df, cluster_col, group1, method, ...) {
    all_clusters <- unique(cluster_df[[cluster_col]])
    other_clusters <- setdiff(all_clusters, group1)
    if (length(other_clusters) == 0) stop("Only one cluster present.")
    cluster_df_temp <- cluster_df
    cluster_df_temp$group_temp <- ifelse(cluster_df_temp[[cluster_col]] == group1, group1, "All_Others")
    res <- run_de_analysis_parallel(
        itx = itx,
        cluster_df = cluster_df_temp,
        cluster_col = "group_temp",
        group1 = group1,
        group2 = "All_Others",
        method = method,
        ...
    )
    res$cluster <- paste0(group1, " vs All_Others")
    return(res)
}