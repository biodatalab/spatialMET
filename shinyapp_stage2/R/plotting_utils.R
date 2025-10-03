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
  need_cols = c('Pixel', 'X', 'Y', 'cluster', 'color')
  if(sum(colnames(pixel_cl) %in% need_cols) == 5){
    pixel_cl = pixel_cl[, need_cols]
  } else{
    stop('Required columns missing in the hcdist output')
  }
  rm(need_cols) # Clean env
  
  # Create pixel index
  # Add initial pixel labels to plot
  # Add leading 'c' to cluster IDs
  pixel_cl[['id']] = c(1:nrow(pixel_cl))
  pixel_cl = setNames(pixel_cl[, colnames(pixel_cl)], nm=tolower(colnames(pixel_cl)))
  pixel_cl[['labels']] = labels
  #pixel_cl[['cluster']] = paste0('c', pixel_cl[['cluster']])
  
  # Extract color palette from hcdist output
  hc_col_pal = pixel_cl[, c('cluster', 'color')]
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
    ggplot2::geom_tile(data=pixel_cl, ggplot2::aes(x=x, y=y, fill=cluster), color=NA, inherit.aes=F, show.legend=F) +
    ggplot2::scale_fill_manual(values=hc_col_pal) +
    ggnewscale::new_scale_fill() +
    ggiraph::geom_tile_interactive(data=pixel_cl, ggplot2::aes(x=x, y=y, data_id=id, fill=labels), color=NA, 
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
#
make_plot_mz = function(mz_dat=NULL, labels=NULL, mz_val=NULL, spot_alpha=NULL){
  
  # Make sure all needed columns are present
  need_cols = c('Pixel', 'X', 'Y', 'mz_intx')
  if(sum(colnames(mz_dat) %in% need_cols) == 4){
    pixel_cl = mz_dat[, need_cols]
  } else{
    stop('Required columns to plot feature intensity are missing.')
  }
  rm(need_cols) # Clean env
  
  # Create pixel index
  # Add initial pixel labels to plot
  # Add leading 'c' to cluster IDs
  mz_dat[['id']] = c(1:nrow(pixel_cl))
  mz_dat = setNames(mz_dat[, colnames(mz_dat)], nm=tolower(colnames(mz_dat)))
  mz_dat[['labels']] = labels
  
  # Color palette for manual annotations
  col_pal = c(Default='#FFFFFF')
  master_col_pal = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  new_labels = unique(mz_dat[['labels']])[ !(unique(mz_dat[['labels']]) %in% names(col_pal)) ]
  if(length(new_labels) > 0){
    new_cols = setNames(master_col_pal[1:length(new_labels)], new_labels)
    col_pal = c(col_pal, new_cols)
  }
  
  gg = ggplot2::ggplot() +
    ggplot2::geom_tile(data=mz_dat, ggplot2::aes(x=x, y=y, fill=mz_intx), color=NA, inherit.aes=F) +
    khroma::scale_fill_sunset(midpoint=median(mz_dat[['mz_intx']]), name='Normalized\nintensity') +
    ggplot2::guides(fill=ggplot2::guide_legend(theme=theme(
      legend.direction="horizontal",
      legend.title.position="top",
      legend.text.position="bottom",
      legend.text=element_text(hjust=1, vjust=0.1, angle=45)
    ))) +
    ggnewscale::new_scale_fill() +
    ggiraph::geom_tile_interactive(data=mz_dat, ggplot2::aes(x=x, y=y, data_id=id, fill=labels), color=NA, 
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

