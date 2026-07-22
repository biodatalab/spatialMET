# R/helpers.R 
# Small utility functions used in multiple places

#' Clean a colour string from various formats
#' @param col character: colour code (e.g., "FF0000", "#FF0000", "xFF0000")
#' @return valid hex colour (#RRGGBB) or #CCCCCC as fallback
clean_color <- function(col) {
    if (is.na(col)) return("#CCCCCC")
    col <- as.character(col)
    col <- sub("^[xX]", "", col)                # remove leading x/X
    if (grepl("^[A-Fa-f0-9]{6}$", col)) col <- paste0("#", col)
    if (grepl("^[A-Fa-f0-9]{3}$", col)) col <- paste0("#", col)
    if (grepl("^#[A-Fa-f0-9]{6}$", col)) return(col)
    return("#CCCCCC")
}

#' Get dominant cluster colors from the cluster data frame
#' @param cluster_vec character vector of cluster labels
#' @param cluster_df data frame with columns hc_manual and color
#' @return named character vector (cluster -> colour)
get_cluster_colors <- function(cluster_vec, cluster_df) {
    unique_clusters <- unique(cluster_vec)
    color_map <- c()
    for (lab in unique_clusters) {
        orig_colors <- cluster_df$color[cluster_df$hc_manual == lab]
        if (length(orig_colors) == 0) {
            color_map[lab] <- "#CCCCCC"
        } else {
            tbl <- table(orig_colors)
            color_map[lab] <- names(tbl)[which.max(tbl)]
        }
    }
    return(color_map)
}

#' Apply intensity transformation (log2 or raw)
#' @param x matrix, DelayedArray, or vector
#' @param method "log2" or "raw"
#' @return transformed object of same class
apply_intensity_transform <- function(x, method = "log2") {
    if (method == "log2") {
        if (inherits(x, "DelayedArray")) {
            x <- log2(x + 1)
        } else {
            x <- log2(x + 1)   # works for matrix and vector
        }
    } else if (method == "raw") {
        # nothing
    } else {
        stop("Unknown transformation method: ", method)
    }
    return(x)
}

#' Default value for NULL or empty
ifnull <- function(x, default) {
    if (is.null(x) || length(x) == 0) default else x
}

#' Safely extract FPCA object from mxfda result
#' @param obj mxfda object
#' @param metric character, e.g. "uni g"
#' @return list with scores, varprop, evalues, or NULL
extract_fpca_safely <- function(obj, metric = "uni g") {
    res <- tryCatch(
        extract_fpca_object(obj, what = paste0(metric, " fpca")),
        error = function(e) NULL
    )
    if (is.null(res)) return(NULL)
    if (!is.list(res)) return(NULL)
    if (is.null(res$scores) || is.null(res$evalues)) return(NULL)
    varprop <- res$evalues / sum(res$evalues)
    list(scores = as.matrix(res$scores),
         varprop = as.numeric(varprop),
         evalues = as.numeric(res$evalues))
}

#' Build a plotly network graph from network_data() output
#' @param dat list with nodes and edges
#' @param show_labels logical: show edge weights?
#' @param label_size numeric: edge label font size
#' @param node_label_size numeric: node label font size
#' @param node_color character
#' @param edge_pos_color character
#' @param edge_neg_color character
#' @return plotly object
plotly_network <- function(dat, show_labels, label_size, node_label_size,
                           node_color, edge_pos_color, edge_neg_color) {
    req(dat)
    
    nodes <- dat$nodes
    nodes$size <- nodes$degree * 2 + 5
    nodes$color <- node_color
    
    edges <- dat$edges
    edges$color <- ifelse(edges$weight > 0, edge_pos_color, edge_neg_color)
    edges$width <- abs(edges$weight) * 3
    
    p <- plotly::plot_ly()
    
    for (i in 1:nrow(edges)) {
        from_id <- edges[i, "from"]
        to_id   <- edges[i, "to"]
        from_x <- nodes$x[nodes$id == from_id]
        from_y <- nodes$y[nodes$id == from_id]
        to_x   <- nodes$x[nodes$id == to_id]
        to_y   <- nodes$y[nodes$id == to_id]
        
        p <- p %>% plotly::add_trace(
            x = c(from_x, to_x, NA),
            y = c(from_y, to_y, NA),
            type = "scatter",
            mode = "lines",
            line = list(color = edges$color[i], width = edges$width[i] / 2),
            hoverinfo = "none",
            showlegend = FALSE
        )
    }
    
    p <- p %>% plotly::add_trace(
        x = nodes$x,
        y = nodes$y,
        type = "scatter",
        mode = "markers+text",
        marker = list(
            size = nodes$size,
            color = nodes$color,
            line = list(color = "black", width = 1)
        ),
        text = nodes$label,
        textposition = "bottom center",
        textfont = list(size = node_label_size),
        hoverinfo = "text",
        hovertext = paste0("Node: ", nodes$label, "<br>Degree: ", nodes$degree)
    )
    
    if (show_labels) {
        mid_anno <- lapply(1:nrow(edges), function(i) {
            from_id <- edges[i, "from"]
            to_id   <- edges[i, "to"]
            mid_x <- mean(nodes$x[nodes$id %in% c(from_id, to_id)])
            mid_y <- mean(nodes$y[nodes$id %in% c(from_id, to_id)])
            list(
                x = mid_x,
                y = mid_y,
                text = as.character(round(edges$weight[i], 3)),
                showarrow = FALSE,
                font = list(size = label_size, color = "black"),
                xanchor = "center",
                yanchor = "bottom"
            )
        })
        p <- p %>% plotly::layout(annotations = mid_anno)
    }
    
    p <- p %>% plotly::layout(
        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        hovermode = "closest",
        title = "Metabolite Correlation Network"
    )
    p
}

#' Compute log2 fold change of a domain vs the rest
#' @param mat matrix (features x pixels)
#' @param labels character vector of domain labels per pixel
#' @param dom character: the domain of interest
#' @return numeric vector of log2 fold changes per feature
compute_log2fc <- function(mat, labels, dom) {
    in_dom  <- labels == dom
    med_in  <- apply(mat[, in_dom,  drop = FALSE], 1, median, na.rm = TRUE)
    med_out <- apply(mat[, !in_dom, drop = FALSE], 1, median, na.rm = TRUE)
    log2((med_in + 1e-6) / (med_out + 1e-6))
}