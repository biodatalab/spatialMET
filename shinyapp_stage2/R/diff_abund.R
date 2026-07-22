##
# @title da_tests: Differential abundance of spatial metabolomics data
#
# @param x a matrix of transformed intensity values per pixel (rows=m/z values and
# columns=pixels)
# @param group1 vector of column indices/names for group 1
# @param group2 vector of column indices/names for group 2
# @param group1_lbl character label for group 1
# @param group2_lbl character label for group 2
# @param method one of "wilcoxon", "ttest", "hellinger"
# @param pval_adj method for p-value adjustment (passed to p.adjust)
#
# @return a data frame with results (metabolite, cluster, avg_log_itx1, avg_log_itx2,
#         log_fc, test_statistic, p_value, adj_p_value)
#
da_tests <- function(x, group1, group2, group1_lbl, group2_lbl, method = "wilcoxon", pval_adj = "fdr") {
    
    # Ensure x is a matrix
    if (!is.matrix(x)) x <- as.matrix(x)
    
    # Convert group indices to numeric positions if they are character names
    if (is.character(group1)) group1 <- which(colnames(x) %in% group1)
    if (is.character(group2)) group2 <- which(colnames(x) %in% group2)
    
    # Pre-compute group means (x is already transformed: raw or log2)
    group1_means <- rowMeans(x[, group1, drop = FALSE], na.rm = TRUE)
    group2_means <- rowMeans(x[, group2, drop = FALSE], na.rm = TRUE)
    log_fc <- group1_means - group2_means  # group1 vs group2 (raw diff or log2 diff)
    
    # Initialize result columns
    test_stat <- rep(NA, nrow(x))
    p_val <- rep(NA, nrow(x))
    
    # Loop through features
    for (i in seq_len(nrow(x))) {
        vals1 <- as.numeric(x[i, group1])
        vals2 <- as.numeric(x[i, group2])
        
        # Remove NA values
        vals1 <- vals1[!is.na(vals1)]
        vals2 <- vals2[!is.na(vals2)]
        
        # Check for sufficient observations (need at least 2 per group for Wilcoxon/t-test)
        if (length(vals1) < 2 || length(vals2) < 2) {
            test_stat[i] <- NA
            p_val[i] <- NA
            next
        }
        
        if (method == "wilcoxon") {
            test <- wilcox.test(vals1, vals2, exact = FALSE)
            test_stat[i] <- test$statistic
            p_val[i] <- test$p.value
        } else if (method == "ttest") {
            test <- t.test(vals1, vals2)
            test_stat[i] <- test$statistic
            p_val[i] <- test$p.value
        } else if (method == "hellinger") {
            # Hellinger distance between two empirical densities
            # Use a simple histogram-based approximation
            all_vals <- c(vals1, vals2)
            breaks <- seq(min(all_vals), max(all_vals), length.out = 50)
            h1 <- hist(vals1, breaks = breaks, plot = FALSE)$density
            h2 <- hist(vals2, breaks = breaks, plot = FALSE)$density
            h1 <- h1 / sum(h1)
            h2 <- h2 / sum(h2)
            hell <- sqrt(0.5 * sum((sqrt(h1) - sqrt(h2))^2))
            test_stat[i] <- hell
            p_val[i] <- NA  # No p-value for Hellinger distance
        }
    }
    
    # Adjust p-values (skip for hellinger where p_val is NA)
    if (method != "hellinger") {
        adj_p <- p.adjust(p_val, method = pval_adj)
    } else {
        adj_p <- rep(NA, nrow(x))
    }
    
    # Build result data frame
    res <- data.frame(
        metabolite = rownames(x),
        cluster = paste(group1_lbl, "vs", group2_lbl),
        avg_log_itx1 = group1_means,
        avg_log_itx2 = group2_means,
        log_fc = log_fc,
        test_statistic = test_stat,
        p_value = p_val,
        adj_p_value = adj_p,
        stringsAsFactors = FALSE
    )
    
    # Arrange by adjusted p-value then absolute logFC (skip if all NA)
    if (!all(is.na(res$adj_p_value))) {
        res <- res %>% dplyr::arrange(adj_p_value, dplyr::desc(abs(log_fc)))
    }
    
    return(res)
}