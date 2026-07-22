##
#' Manual ROI selection tool via shiny
#' ADAPTED FROM: STUtility
#' This function takes a seurat object with stored image locations and opens up the manual selection tool in the default browser
#'
#' @param coords Data frame with coordinates. Three columns: pixel name, y and x
#' @param intx Intensity matrix
#' @param mz The m/z value (metabolite) to plot in the background
#'
#' @export
#'
manual_annotation_2 = function(coords=NULL, intx=NULL, mz=NULL){
  
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
                 #selectInput(inputId = "sampleInput", label = "Select sample", choices = sampleChoice, selected = "1"),
                 #shiny::hr(),
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
      #ggiraph::ggiraphOutput("Plot1", width="100%", height=paste0(res, "px"))
      ggiraph::girafeOutput("Plot1", width="100%", height=paste0(res, "px"))
    )
  )

  # ===================== Server ================================
  server = function(input, output, session){

    df_tmp = reactiveValues(label=rep('Default', nrow(coords)),
                            id=c(1:nrow(coords)),
                            sample=rep('1', nrow(coords)))

    if(!is.null(intx) & mz_plot %in% intx[[1]]){
      intx_mz = as.data.frame(t(intx[ intx[[1]] == mz_plot, -1]))
      colnames(intx_mz) = 'itx_val'
      intx_mz[['pixname']] = rownames(intx_mz)
      intx_mz = intx_mz[, c('pixname', 'itx_val')]
    } else{
      intx_mz = NULL
    }
    
    # rv <- reactiveValues(ann=NULL)
    # observeEvent("1", {
    #   rv$ann <- create_annotation(coords, intx_mz)
    # })
    
    output$Plot1 <- ggiraph::renderGirafe({
      
      withProgress(message="Updating plot", value=0, {
        x_ggr = ggiraph::girafe(ggobj=make_plot(res=res,
                                                spotAlpha=input$alphaValue,
                                                Labels=df_tmp$label,
                                                SpotSize=input$spotSize,
                                                intx_mz=intx_mz,
                                                coordinates=coords
                                                #ann = rv$ann
                                                #spotAlpha=0.1, ###TEST
                                                #Labels=rep('Default', nrow(coords)), ###TEST
                                                #SpotSize=1, ###TEST
                                                #ann=create_annotation(coords, intx_mz) ###TEST
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
        stopApp(returnValue=coords)
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
#' @param spotAlpha geom_point opacity
#' @param Labels labels of the spots (Usually "Default")
#' @param SpotSize geom_point size
#' @param res resolution of the image
#' @param ann Annotation object (contains intensity values, coordinates, and other parameters)
#'
#' @importFrom magick image_read image_scale image_info geometry_size_pixels
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous theme element_rect element_blank theme_minimal
#'
#' @keywords internal

make_plot = function (spotAlpha, Labels, SpotSize, res, intx_mz, coordinates){
  
  px.ids = colnames(coordinates)
  coordinates[['id']] = c(1:nrow(coordinates))
  coordinates <- setNames(coordinates[, c(px.ids, "id")], nm=c("pixname", "y", "x", "id"))
  coordinates[['Labels']] = Labels
  
  coords_match_key = colnames(coordinates)[1]
  coordinates = intx_mz %>% dplyr::left_join(., coordinates, by='pixname')
  
  mid = median(intx_mz[['itx_val']])
  
  col_pal = c("gray90", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  gg = ggplot() +
    #ann +
    ggiraph::geom_tile_interactive(data=coordinates, aes(x=x, y=y, data_id=id, color=Labels, fill=itx_val), 
                                   #shape=0, 
                                   stroke=SpotSize, alpha=spotAlpha, inherit.aes=F) +
    coord_fixed() +
    scale_color_discrete(type=col_pal) +
    scale_fill_gradient2(midpoint=mid, low="blue", mid="white", high="red") +
    #scale_fill_gradient(low="blue", high="red") +
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
#' @param coords Coordinates data frame
#' @param intx_mz A two column data frame with pixel names matching coords and the itensity values
#' 
#' @importFrom stats setNames
#' 
create_annotation = function(coords=NULL, intx_mz=NULL){
  
  px.ids = colnames(coords)
  coordinates = coords
  coordinates[['id']] = c(1:nrow(coordinates))
  coordinates <- setNames(coordinates[, c(px.ids, "id")], nm=c("pixname", "y", "x", "id"))

  if(!is.null(intx_mz)){
    coords_match_key = colnames(coordinates)[1]
    intx_mz = intx_mz %>% dplyr::left_join(., coordinates, by='pixname')
    #mz_layer = ggplot2::geom_point(data=intx_mz, aes(x=x, y=y, fill=itx_val), stroke=NA, shape=22, inherit.aes=F)
    mz_layer = ggplot2::geom_raster(data=intx_mz, aes(x=x, y=y, fill=itx_val), inherit.aes=F)
    median_itx = median(intx_mz[['itx_val']])
  } else {
    mz_layer <- NULL
    median_itx = 0
  }

  return(list(coordinates, median_itx, mz_layer))
}

