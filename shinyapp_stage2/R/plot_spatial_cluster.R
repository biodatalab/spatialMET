##
#' @title plot_spatial_clusters: Plot cluster memberships of ST spots
#' @description Generates a plot of the location of spots within an spatial array,
#' and colors them according to spatially-weighted hierarchical clustering assignments.
#' @details
#' The function takes an STList with cluster memberships and plots the spots with
#' colors indicating the cluster they belong to. Optionally, the user can annotate
#' tumor/stroma compartments if ESTIMATE scores are available.
#'
#' @param x an STList with hierarchical cluster memberships.
#' @param coords a data frame with three columns. The first column containing the names of
#' the pixels matching column names of the intensity matrox (`x`). The second and
#' third columns containing the x and y coordinates of the pixels
#' @param cluster_df a data frame. The output of `spatial_clust`
#' @param ks the k values to plot or `dtc` for dynamictreecluster
#' @param ws the spatial weights to plot
#' @param plot_meta a column name from the output of `spatial_clust`
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @param ptsize a number specifying the size of the points. Passed to `size` aesthetic.
#' @param txsize a number specifying the size of the font. Passed to `size` argument from `element_text`.
#' @return a list with the requested plots.
#'
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#
#
plot_spatial_clusters = function(coords=NULL, cluster_df=NULL, ks='dtc', ws=NULL, 
                                 deepSplit=NULL, plot_meta=NULL, color_pal=NULL, 
                                 ptsize=NULL, txsize=NULL, cluster_filter=NULL) {
    
    # Check coordinates and cluster assignments were provided
    if(is.null(coords) | is.null(cluster_df)){
        stop('Missing coordinates and/or cluster assignments.')
    }
    
    # Define columns to plot
    if(is.null(plot_meta)){
        plot_meta = grep('^spatial_clust_spw', colnames(cluster_df), value=T)
        
        if(!is.null(ws)){
            if(any(ws == 0)){
                ws_tmp = ws[ws != 0]
                plot_meta_tmp = grep('spatial_clust_spw0_|spatial_clust_spw0$', plot_meta, value=T)
                if(length(ws_tmp) > 0){
                    plot_meta = c(plot_meta_tmp, grep(paste0('spatial_clust_spw', ws_tmp, collapse='|'), plot_meta, value=T))
                } else {
                    plot_meta = plot_meta_tmp
                }
                rm(ws_tmp, plot_meta_tmp)
            } else {
                plot_meta = grep(paste0('spatial_clust_spw', ws, collapse='|'), plot_meta, value=T)
            }
        }
        
        if(ks[1] != 'dtc'){
            plot_meta = grep(paste0('_k', ks,'$', collapse='|'), plot_meta, value=T)
        } else if(ks[1] == 'dtc'){
            if(!is.null(deepSplit)){
                plot_meta = grep(paste0('_dspl', stringr::str_to_title(as.character(deepSplit)), '$', collapse='|'), plot_meta, value=T)
            } else {
                plot_meta = grep('_dspl', plot_meta, value=T)
            }
        } else {
            stop('Specify one or several k values to plot, or use ks=\'dtc\' (default).')
        }
    }
    
    # Check that the meta data column exists
    if(length(grep(paste0('^', plot_meta, '$', collapse='|'),  colnames(cluster_df))) == 0){
        stop('No metadata column or clustering parameters were specified. Or specified parameters do not exist in metadata.')
    }
    
    plot_list = list()
    
    # Define size of points
    if(is.null(ptsize)){
        ptsize = 0.5
    }
    
    # Define size of text
    if(is.null(txsize)){
        txsize = 12
    }
    
    for(metacol in plot_meta){
        # Set default color if NULL input
        if(is.null(color_pal)){
            color_pal = 'light'
            if(is.numeric(cluster_df)){
                color_pal = 'sunset'
            }
        }
        
        # Join coordinates and cluster data
        df_tmp = cluster_df %>%
            dplyr::left_join(coords %>% dplyr::rename(pixname=1, ypos=2, xpos=3), ., by='pixname') %>%
            dplyr::select(pixname, xpos, ypos, meta:=!!metacol)
        
        # ------------------------------------------------------------
        # Apply cluster filter if provided
        # ------------------------------------------------------------
        if(!is.null(cluster_filter)){
            # Convert filter to character for safe comparison
            filter_vals <- as.character(cluster_filter)
            df_tmp <- df_tmp %>% dplyr::filter(as.character(meta) %in% filter_vals)
            # If no rows remain, skip this metacol
            if(nrow(df_tmp) == 0){
                warning(paste("No pixels found for cluster_filter =", 
                              paste(cluster_filter, collapse=", "), "in column", metacol))
                next
            }
        }
        
        if(!is.numeric(df_tmp[['meta']]) & length(unique(df_tmp[['meta']])) < 100){
            # Convert meta data to factor
            df_tmp = df_tmp %>%
                dplyr::mutate(meta=tidyr::replace_na(as.character(meta), 'No_Data')) %>%
                dplyr::mutate(meta=as.factor(meta))
            
            # Create color palette.
            meta_cols = color_parse(color_pal, n_cats=length(unique(df_tmp[['meta']])))
            names(meta_cols) = unique(df_tmp[['meta']])
            if(any(names(meta_cols) == 'No_Data')){
                meta_cols[names(meta_cols) == 'No_Data'] = 'gray50'
            }
        }
        
        # Prepare titles for plots
        if(grepl('^spatial_clust_', metacol)){
            title_w = as.character(stringr::str_extract(metacol, paste0("spw0\\.?[0-9]*"))) %>% gsub('spw', '', .)
            if(grepl('_k[0-9]+$', metacol)){
                title_k = as.character(stringr::str_extract(metacol, paste0("_k[0-9]+"))) %>% gsub('_k', '', .)
                title_p = paste0("k=", title_k, "\nspatial weight=", title_w)
            } else if(grepl('_dspl[\\.0-9TrueFalse]+$', metacol)){
                title_dspl = stringr::str_extract(metacol, '[\\.0-9]+$|True$|False$')
                title_p = paste0("DynamicTreeCut; deepSplit=", title_dspl, "\nSpatial weight=", title_w)
            }
            title_leg = 'Clusters'
        } else {
            title_leg = metacol
        }
        
        # Append filter info to title if applied
        if(!is.null(cluster_filter)){
            title_p = paste0(title_p, "\n(Showing only cluster(s): ", 
                             paste(cluster_filter, collapse=", "), ")")
        }
        
        # Create plot
        p = ggplot2::ggplot(df_tmp) +
            ggplot2::geom_point(ggplot2::aes(x=xpos, y=ypos, color=meta), size=ptsize)
        
        # Assign color palette
        if(is.factor(df_tmp[['meta']])){
            p = p + ggplot2::scale_color_manual(values=c(meta_cols))
        } else {
            if(!is.null(color_pal)){
                khroma_cols = khroma::info()
                khroma_cols = khroma_cols$palette
                if(color_pal %in% khroma_cols){
                    p = p + khroma::scale_color_picker(palette=color_pal)
                } else {
                    p = p + ggplot2::scale_color_distiller(palette=color_pal)
                }
            }
        }
        
        if(!is.numeric(df_tmp[['meta']])){
            p = p +
                ggplot2::guides(color=guide_legend(override.aes=list(size=ptsize+1)))
        }
        
        p = p +
            labs(color=title_leg, title=title_p) +
            ggplot2::coord_fixed(ratio=1) +
            ggplot2::theme_void() +
            ggplot2::theme(legend.title=element_text(size=txsize),
                           plot.title=element_text(size=txsize+2))
        
        plot_list[[metacol]] = p
    }
    
    return(plot_list)
}