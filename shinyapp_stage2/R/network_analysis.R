# R/network_analysis.R
library(igraph)
library(visNetwork)
library(tidygraph)
library(ggraph)

#' Build Metabolite Correlation Network
#' 
#' @param itx Intensity matrix (features x samples/pixels)
#' @param top_n Top variable features to include
#' @param corr_method "pearson" or "spearman"
#' @param corr_threshold Absolute correlation threshold
#' @param min_degree Minimum degree for nodes to keep
#' Build Metabolite Correlation Network
build_metabolite_network <- function(itx, 
                                     pixel_clusters = NULL,
                                     selected_domain = NULL,     # NEW: single domain
                                     top_n = 300,
                                     top_n_per_domain = 50,
                                     use_per_domain = FALSE,
                                     corr_method = "pearson",
                                     corr_threshold = 0.6,
                                     min_degree = 2) {
    
    if (is.null(itx) || nrow(itx) == 0) {
        return(list(nodes = data.frame(), edges = data.frame(), message = "No data"))
    }
    
    # === DOMAIN-SPECIFIC FEATURE SELECTION ===
    if (use_per_domain && !is.null(pixel_clusters) && !is.null(selected_domain)) {
        
        pix_in_dom <- pixel_clusters$pixel_id[pixel_clusters$hc_manual == selected_domain]
        pix_in_dom <- intersect(pix_in_dom, colnames(itx))
        
        if (length(pix_in_dom) == 0) {
            return(list(nodes = data.frame(), edges = data.frame(), 
                        message = paste("No pixels found for domain:", selected_domain)))
        }
        
        sub_itx <- itx[, pix_in_dom, drop = FALSE]
        vars <- apply(sub_itx, 1, sd, na.rm = TRUE)
        top_idx <- order(vars, decreasing = TRUE)[1:min(top_n_per_domain, length(vars))]
        selected_mz <- rownames(itx)[top_idx]
        
        itx <- itx[selected_mz, , drop = FALSE]
        message_text <- paste("Domain:", selected_domain, "| Top", length(selected_mz), "features")
        
    } else {
        # Global top variable
        if (nrow(itx) > top_n) {
            vars <- apply(itx, 1, sd, na.rm = TRUE)
            keep <- order(vars, decreasing = TRUE)[1:top_n]
            itx <- itx[keep, , drop = FALSE]
        }
        message_text <- paste("Global top", nrow(itx), "features")
    }
    
    # === Build Correlation Network ===
    cor_mat <- cor(t(itx), method = corr_method, use = "pairwise.complete.obs")
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 0
    
    edge_list <- which(abs(cor_mat) >= corr_threshold, arr.ind = TRUE)
    if (nrow(edge_list) == 0) {
        return(list(nodes = data.frame(), edges = data.frame(), 
                    message = "No edges above correlation threshold"))
    }
    
    edges <- data.frame(
        from = rownames(cor_mat)[edge_list[,1]],
        to = colnames(cor_mat)[edge_list[,2]],
        weight = cor_mat[edge_list],
        stringsAsFactors = FALSE
    )
    
    g <- graph_from_data_frame(edges, directed = FALSE)
    V(g)$degree <- degree(g)
    V(g)$strength <- strength(g, weights = abs(E(g)$weight))
    
    g <- induced_subgraph(g, V(g)$degree >= min_degree)
    
    if (vcount(g) == 0) {
        return(list(nodes = data.frame(), edges = data.frame(), 
                    message = "No nodes after degree filtering"))
    }
    
    nodes <- data.frame(
        id = V(g)$name,
        label = V(g)$name,
        value = V(g)$strength,
        title = paste0(V(g)$name, "<br>Degree: ", V(g)$degree),
        stringsAsFactors = FALSE
    )
    
    edges_df <- igraph::as_data_frame(g)
    edges_df$width <- abs(edges_df$weight) * 3
    
    list(
        nodes = nodes,
        edges = edges_df,
        graph = g,
        message = paste(message_text, "| Network:", vcount(g), "nodes,", ecount(g), "edges")
    )
}