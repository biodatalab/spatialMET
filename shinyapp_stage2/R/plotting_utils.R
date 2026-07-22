##
# Plotting functions
#

##
# make_plot_hcdist: Used for the manual annotation tool, returns spatial plot colored by hcdist cluster
#
# @param pixel_cl
# @param labels labels of the spots (Usually "Default")
# @param spot_alpha geom_point opacity
#
make_plot_hcdist = function(pixel_cl=NULL, labels=NULL, spot_alpha=NULL){

  # Convert pixel_cl to R data frame if data.table
  if(any(class(pixel_cl) == 'data.table')){
    pixel_cl = as.data.frame(pixel_cl)
  }

  # Make sure all needed columns are present
  need_cols = c('row_id', 'pixel_id', 'x_coord', 'y_coord', 'hc_orig', 'color')
  if(sum(colnames(pixel_cl) %in% need_cols) == 6){
    pixel_cl = pixel_cl[, need_cols]
  } else{
    stop('Required columns missing in the hcdist output') # NEED TO ADD CASE FOR NON-HCDIST INPUT CASE
  }
  rm(need_cols) # Clean env

  # Add initial manual pixel labels to plot
  pixel_cl[['labels']] = labels

  # Extract color palette from hcdist output
  hc_col_pal = pixel_cl[, c('hc_orig', 'color')]
  hc_col_pal = unique(hc_col_pal)
  hc_col_pal[['color']] = gsub('^x', '#', hc_col_pal[['color']])
  hc_col_pal = setNames(hc_col_pal[[2]], hc_col_pal[[1]])

  # Color palette for manual annotations
  col_pal = hc_col_pal
  master_col_pal = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  new_labels = unique(pixel_cl[['labels']])[ !(unique(pixel_cl[['labels']]) %in% names(hc_col_pal)) ]
  if(length(new_labels) > 0){
    new_cols = setNames(master_col_pal[1:length(new_labels)], new_labels)
    col_pal = c(hc_col_pal, new_cols)
  }

  gg = ggplot2::ggplot() +
    ggplot2::geom_tile(data=pixel_cl, ggplot2::aes(x=x_coord, y=y_coord, fill=hc_orig), color=NA, inherit.aes=F, show.legend=F) +
    ggplot2::scale_fill_manual(values=hc_col_pal) +
    ggnewscale::new_scale_fill() +
    ggiraph::geom_tile_interactive(data=pixel_cl, ggplot2::aes(x=x_coord, y=y_coord, data_id=row_id, fill=labels), color=NA,
                                   alpha=spot_alpha, inherit.aes=F) +
    ggplot2::scale_fill_manual(values=col_pal, name='Manual\nannotation') +
    ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(alpha=1))) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="black", colour=NA),
                   legend.position="bottom",
                   axis.text=ggplot2::element_blank(),
                   axis.title=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   panel.grid=ggplot2::element_blank())

  return(gg)
}


##
# make_plot_mz: Used for the manual annotation tool, returns spatial plot showing intensities for a given m/z
#
# @param mz_dat
# @param mz_val name of feature
# @param spot_alpha geom_point opacity
make_plot_mz_optimized <- function(
        mz_dat,
        mz_val = NULL,
        spot_alpha = 0.7,
        spot_color = "#1f78b4",
        contrast = 1
) {
    
    if (!data.table::is.data.table(mz_dat)) {
        mz_dat <- data.table::as.data.table(mz_dat)
    }
    
    # Ensure labels column exists
    if (!"labels" %in% names(mz_dat)) {
        mz_dat[, labels := "Default"]
    }
    
    # Intensity gradient: white → spot_color
    fill_scale <- scale_fill_gradient(
        low = "white", 
        high = spot_color,
        name = "Intensity",
        limits = c(0, 1),
        oob = scales::squish
    )
    
    # Annotation color palette (single color)
    unique_labels <- unique(mz_dat$labels)
    ann_pal <- setNames(rep(spot_color, length(unique_labels)), unique_labels)
    
    # Build plot
    p <- ggplot() +
        # === Intensity layer (fast raster) ===
        geom_raster(
            data = mz_dat,
            aes(x = x_coord, y = y_coord, fill = mz_intx)
        ) +
        fill_scale +
        
        # === Switch to new fill scale for annotations ===
        ggnewscale::new_scale_fill() +
        
        # === Interactive annotation overlay ===
        ggiraph::geom_tile_interactive(
            data = mz_dat,
            aes(
                x = x_coord, 
                y = y_coord,
                data_id = row_id,
                fill = labels
            ),
            color = NA,
            alpha = spot_alpha,
            inherit.aes = FALSE
        ) +
        
        scale_fill_manual(
            values = ann_pal,
            name = "Manual annotation",
            guide = guide_legend(override.aes = list(alpha = 1))
        ) +
        
        coord_fixed() +
        scale_y_reverse() +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "black", colour = NA),
            legend.position = "bottom",
            legend.box = "vertical",
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()
        )
    
    p
}