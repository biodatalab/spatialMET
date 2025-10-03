##
#' Manual ROI selection tool via shiny
#' ADAPTED FROM: STUtility
#' This function takes outputs from hcdist (Eric Welsh) and opens up the manual selection tool in the default browser
#'
#' @param pixel_clusters the "pixel_clusters.txt" file resulting from hcdist
#' @param intx Intensity matrix used as input for hcdist
#' @param mz The m/z value (metabolite) to plot in the background
#'
#' @export
#'
manual_annotation_hcdist = function(pixel_clusters=NULL, intx=NULL, mz_plot=NULL){
  
  library('shiny')
  
  # Input image type, e.g. "masked" or "raw". Keep NULL for now
  type=NULL
  # Keep option for future implementation with tissue image (width in pixels)
  res=1000
  
  # ===================== UI =======================
  ui = pageWithSidebar(
    headerPanel("Manual selection"),
    
    sidebarPanel(width=3,
                 actionButton(inputId="info", label="Instructions"),  #icon = shiny::icon("info", lib="glyphicon")),
                 shiny::hr(),
                 textInput(inputId="labelInput", label="Choose ROI/annotation name", value="", placeholder="Default"),
                 shiny::hr(),
                 sliderInput(inputId="alphaValue", label="Opacity [0-1]", min=0, max=1, value=0.2, step=0.1),
                 shiny::hr(),
                 sliderInput(inputId="spotSize", label="Pixel sixe [0.1-1]", min=0, max=1, value=0.5, step=0.1),
                 shiny::hr(),
                 actionButton(inputId="confirm", label="Confirm ROI selection"),
                 shiny::hr(),
                 actionButton(inputId="stopApp", label="Quit annotation tool")
    ),
    mainPanel(
      ggiraph::girafeOutput("Plot1", width="100%", height=paste0(res, "px"))
    )
  )
  
  # ===================== Server ================================
  server = function(input, output, session){
    
    df_tmp = reactiveValues(label=rep('Default', nrow(pixel_clusters)),
                            id=c(1:nrow(pixel_clusters)),
                            sample=rep('1', nrow(pixel_clusters)))
    
    if(!is.null(intx) & mz_plot %in% intx[[1]]){
      intx_mz = as.data.frame(t(intx[ intx[[1]] == mz_plot, -1]))
      colnames(intx_mz) = 'itx_val'
      intx_mz[['pixel']] = rownames(intx_mz)
      intx_mz = intx_mz[, c('pixel', 'itx_val')]
    } else{
      intx_mz = NULL
    }
    
    output$Plot1 <- ggiraph::renderGirafe({
      
      withProgress(message="Updating plot", value=0, {
        x_ggr = ggiraph::girafe(ggobj=make_plot(spot_alpha=input$alphaValue,
                                                labels=df_tmp$label,
                                                pixel_cl=pixel_clusters
        ), width_svg=12, height_svg=10)
        
        x_ggr = ggiraph::girafe_options(x_ggr,
                                        ggiraph::opts_zoom(max = 5),
                                        ggiraph::opts_selection(type = "multiple",
                                                                css = "fill:cyan;stroke:black;opacity:0.7;"))
        return(x_ggr)
      })
    })
    
    observeEvent(input$confirm, {
      ids.selected <- as.numeric(input$Plot1_selected)
      df_tmp$label[which(df_tmp$id %in% ids.selected)] = input$labelInput
      session$sendCustomMessage(type='Plot1_set', message = character(0))
    })
    
    observe({
      if(input$stopApp > 0){
        print("Stopped")
        coords[['ROI_labels']] = df_tmp[['label']]
        stopApp(returnValue=pixel_clusters)
      }
    })
    
    observeEvent(input$info, {
      showModal(modalDialog(
        title = "Instructions",
        HTML("1. Specifiy label you want to use<br>",
             "2. Use the blue(select) lasso to label the caputure-spots<br>",
             "3. Press Confirm to set the labels<br>",
             "4. Repeat 1-4 until all labels are set<br>",
             "5. Close the shiny tool to return"),
        easyClose = TRUE,
        footer = NULL
      ))
    })
    
  } # CLOSE SERVER
  
  runApp(list(ui=ui, server=server), launch.browser=T)
}

#' Used for the manual annotation tool, returns plot and coordinate IDs for selected sample
#'
#' @param pixel_cl
#' @param spot_alpha geom_point opacity
#' @param labels labels of the spots (Usually "Default")
#'
#' @importFrom magick image_read image_scale image_info geometry_size_pixels
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous theme element_rect element_blank theme_minimal
#'
#' @keywords internal
#
make_plot = function(pixel_cl, spot_alpha, labels){
  
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
  pixel_cl[['cluster']] = paste0('c', pixel_cl[['cluster']])
  
  # Extract color palette from hcdist output
  hc_col_pal = pixel_cl[, c('cluster', 'color')]
  hc_col_pal = unique(hc_col_pal)
  hc_col_pal[['color']] = gsub('^x', '#', hc_col_pal[['color']])
  hc_col_pal = setNames(hc_col_pal[[2]], hc_col_pal[[1]])
  
  # Color palette for manual annotations
  col_pal = c("gray90", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  gg = ggplot() +
    ggiraph::geom_tile_interactive(data=pixel_cl, aes(x=x, y=y, data_id=id, color=labels, fill=cluster), 
                                   alpha=spot_alpha, inherit.aes=F) +
    coord_fixed() +
    scale_color_discrete(type=col_pal) +
    scale_fill_manual(values=hc_col_pal) +
    labs(color='ROIs', fill='hcdist') +
    theme_minimal() +
    theme(panel.background=element_rect(fill="transparent", colour=NA),
          legend.position="bottom",
          axis.text=element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank())
  
  return(gg)
}

#' Create an annotation
#' 
#' @param pixel_cl
#' @param intx_mz A two column data frame with pixel names matching coords and the intensity values
#' 
#' @importFrom stats setNames
#' 
create_annotation = function(pixel_cl=NULL, intx_mz=NULL){
  
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
  pixel_cl[['id']] = c(1:nrow(pixel_cl))
  pixel_cl = setNames(pixel_cl[, colnames(pixel_cl)], nm=tolower(colnames(pixel_cl)))
  pixel_cl = pixel_cl[, c('pixel', 'x', 'y', 'id')]
  
  if(!is.null(intx_mz)){
    intx_mz = intx_mz %>% dplyr::left_join(., pixel_cl, by='pixel')
    mz_layer = ggplot2::geom_raster(data=intx_mz, aes(x=x, y=y, fill=itx_val), inherit.aes=F)
    median_itx = median(intx_mz[['itx_val']])
  } else {
    mz_layer <- NULL
    median_itx = 0
  }
  
  return(list(pixel_cl, median_itx, mz_layer))
}

