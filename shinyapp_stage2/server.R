##
# Shiny app for downstream analysis of hcdist results (Eric Welsh's algorithm)
#

# Allow large file upload
options(shiny.maxRequestSize=10000*1024^2)

# ================================= Server =====================================
shinyServer(function(input, output, session){

  # Load hcdist input file from user
  hcdist_res = eventReactive(input$hcdist_user_file, {
    req(valid_hcdist_file())
    infile = input$hcdist_user_file

    # Read cluster file
    cat('READING CLUSTER TABLE\n')
    pixel_clusters = data.table::fread(infile$datapath,
                                       check.names=FALSE, data.table=FALSE,
                                       colClasses=list(character=c("cluster")))

    # Define column names expected if hcdist clusters were input
    hcdist_headers = c("Pixel", "X", "Y", "tree order", "cluster", "cluster_size", "color")

    # Check if file contains hcdist's headers. If not process as generic
    ## X and Y coordinates need to be integers (geom_tile)!!!
    if(sum(colnames(pixel_clusters) %in% hcdist_headers) != 7){
      pixel_clusters = pixel_clusters[, c(1:4)]
      colnames(pixel_clusters) = c("pixel_id", "x_coord", "y_coord", "hc_orig")

      # Create color palette and add to dataframe
      set.seed(12345)
      col_pal = sample(khroma::color("smoothrainbow", force=TRUE)(length(unique(pixel_clusters[['hc_orig']]))))
      col_pal = data.frame(cluster=unique(pixel_clusters[['hc_orig']]), color=col_pal)
      pixel_clusters = pixel_clusters %>% dplyr::left_join(., col_pal, by='hc_orig')

    } else if(sum(colnames(pixel_clusters) %in% hcdist_headers) == 7){
      # Rename columns if resulting from hcdist
      pixel_clusters = pixel_clusters %>%
        dplyr::select(c("pixel_id"="Pixel",
                        'x_coord'="X", 'y_coord'="Y",
                        'hc_orig'="cluster", "color"))
    }

    # Force clusters as character and coordinates as integers
    pixel_clusters[['hc_orig']] = as.character(pixel_clusters[['hc_orig']])
    pixel_clusters[['x_coord']] = as.integer(pixel_clusters[['x_coord']])
    pixel_clusters[['y_coord']] = as.integer(pixel_clusters[['y_coord']])

    # Add columns for manual ROI annotations
    pixel_clusters[['hc_manual']] = pixel_clusters[['hc_orig']]
    pixel_clusters[['mz_manual']] = rep('Default', nrow(pixel_clusters))

    # Sort coordinates to help in pixel/neighborhoods collapsing
    # Add row IDs for manual ROI selection
    pixel_clusters = pixel_clusters %>%
      dplyr::arrange(x_coord, y_coord) %>%
      tibble::add_column(row_id=as.integer(c(1:nrow(.))), .before=1)

    return(pixel_clusters)
  })

  # Prepare reactive value to store data frame with HCDIST and manual annotations
  pixel_clusters = reactiveVal()
  observe({
    pixel_clusters(hcdist_res())
  })

  # If too many pixels, distance matrix eats all memory.
  ## So, find neighborhoods to collapse data
  pixel_collapse = reactive({
    req(valid_hcdist_file())
    coords_df = pixel_clusters()

    min_nb = as.integer(input$collapse_par)
    if(nrow(coords_df) >= 50000 | !is.integer(min_nb) | length(min_nb) == 0){
      min_nb = 21
    }
assign('min_nb', min_nb, envir = .GlobalEnv)
assign('coords_df', coords_df, envir = .GlobalEnv)
    nbs_ls = NULL
    if(min_nb != ""){
      cat('COLLAPSING PIXELS\n')
      withProgress(message="Collapsing pixels to reduce memory usage...", value=0, {
        coords_df = coords_df[, c('pixel_id', 'x_coord', 'y_coord')]
        nbs_ls = find_adj_neighbors(coords_df=coords_df)
        collapse_itx = itx_res()
        collapse_itx = lapply(nbs_ls, function(i){
          itx_tmp = apply(collapse_itx[, i, drop=FALSE], 1, median)
        })
      })
    }
 assign('nbs_ls', nbs_ls, envir = .GlobalEnv)
    return(nbs_ls)
  })

  neighbors = reactiveVal()
  observe({
    neighbors = pixel_collapse()
assign('neighbors', neighbors, envir = .GlobalEnv)
  })

  # Load intensity matrix file from user
  itx_res = reactive({
    req(valid_itx_file())

    withProgress(message="Reading intensity data...", value=0, {
      infile = input$itx_user_file
      if(is.null(infile)){
        intx = NULL
      } else{
        cat('READING INTENSITY TABLE\n')
        intx = vroom::vroom(infile$datapath, guess_max=10, show_col_types=FALSE)
        intx = tibble::column_to_rownames(.data=intx, var=colnames(intx)[1])
        # Convert the data frame to a DelayedArray object
        intx = DelayedArray::DelayedArray(as.matrix(intx))
      }
    })

    return(intx)
  })

  # Extract top variable metabolites
  top_var = reactive({
    itx_tmp = itx_res()
    if(!is.null(itx_tmp)){ # To avoid error when loading selectizeInput before data (at app start)
      var_mz = sort(apply(itx_tmp, 1, sd), decreasing=T)
      return(names(var_mz))
    }
  })

  # Element providing the list of metabolites to select from
  observe({
    updateSelectizeInput(session, 'mz_selected', choices=top_var(), server=T)
  })

  # Element providing the list of clusters to select cluster 1 for DA
  observe({
    if(input_validate()){
      cluster_opts = hcdist_res()
      cluster_opts = unique(cluster_opts[[input$annotation_test]])
      updateSelectizeInput(session, 'group1_selected', choices=cluster_opts, server=T, selected=cluster_opts[1])
    }
  })

  # Element providing the list of clusters to select cluster 2 for DA
  observe({
    if(input_validate()){
      cluster_opts = hcdist_res()
      cluster_opts = unique(cluster_opts[[input$annotation_test]])
      updateSelectizeInput(session, 'group2_selected', choices=cluster_opts, server=T, selected=cluster_opts[2])
    }
  })

  # Element providing the list of reference clusters to spatial gradients
  observe({
    if(input_validate()){
      cluster_opts = hcdist_res()
      cluster_opts = unique(cluster_opts[[input$annotation_test_gradients]])
      updateSelectizeInput(session, 'ref_group', choices=cluster_opts, server=T, selected=cluster_opts[3])
    }
  })


  ############################### SPATIAL DOMAINS ##############################
  # Generate plot showing clusters/domains
  base_hc_plot = reactive({
    dat_tmp = pixel_clusters()
    labs = dat_tmp[['hc_manual']]
    alpha_val = as.numeric(input$alpha_value_hc)
    hc_p = make_plot_hcdist(pixel_cl=dat_tmp,
                            labels=labs,
                            spot_alpha=alpha_val)

    return(hc_p)
  })

  # Generate plot of hcdist clusters
  output$hc_plot = ggiraph::renderGirafe({
    withProgress(message="Updating plot...", value=0, {
      x_ggr = ggiraph::girafe(ggobj=base_hc_plot(),
                              width_svg=12, height_svg=12,
                              ggiraph::opts_selection(type='single'))

      x_ggr = ggiraph::girafe_options(x_ggr, ggiraph::opts_zoom(max=5),
                                      ggiraph::opts_selection(type="multiple",
                                                              css="fill:cyan;stroke:black;opacity:0.7;"))

      return(x_ggr)
    })
  })

  # Generate plot to present tissue image (hcdist module)
  output$img_hc_plot = renderImage({
    if(!is.null(input$img_user_file)){
      list(src=input$img_user_file$datapath)
    }
  }, deleteFile=F)

  # Spatial domain plot PDF download
  output$download_hc_pdf = downloadHandler(
    filename = function(){
      paste0("spatial_domains_", Sys.Date(), ".pdf")
    },
    content = function(file){
      ggsave(file, plot=base_hc_plot(), device="pdf", width=8, height=6)
    }
  )

  # Spatial domain plot PNG download
  output$download_hc_png = downloadHandler(
    filename = function(){
      paste0("spatial_domains_", Sys.Date(), ".png", sep="")
    },
    content = function(file){
      ggsave(file, plot=base_hc_plot(), device="png", width=8, height=6, dpi=300)
    }
  )

  # Download annotations from hcdist page
  output$export_df_hc = downloadHandler(
    filename = function(){
      paste0("manual_annotations_", Sys.Date(), ".csv")
    },
    content = function(file){
      df_tmp = pixel_clusters()
      df_tmp = data.frame(id=df_tmp[['pixel_id']],
                          x_coord=df_tmp[['x_coord']],
                          y_coord=df_tmp[['y_coord']],
                          manual_from_hcdist=df_tmp[['hc_manual']],
                          manual_from_intx=df_tmp[['mz_manual']])

      write.csv(df_tmp, file, quote=T, row.names=F)
    }
  )


  ################################# INTENSITIES ################################
  # Generate plot of mz values
  base_mz_plot = reactive({
    # Feature to plot
    mz_val = input$mz_selected
    # Prepare data frame
    mtx_tmp = as.data.frame(t(itx_res()[mz_val, , drop=FALSE]))
    colnames(mtx_tmp) = 'mz_intx'
    mtx_tmp[['pixel_id']] = rownames(mtx_tmp)

    dat_tmp = pixel_clusters()
    labs = dat_tmp[['mz_manual']]
    alpha_val = as.numeric(input$alpha_value_mz)

    mtx_tmp = mtx_tmp %>% dplyr::left_join(dat_tmp, ., by='pixel_id')
    mz_p = make_plot_mz(mz_dat=mtx_tmp,
                        labels=labs,
                        mz_val=mz_val,
                        spot_alpha=alpha_val)

    return(mz_p)
  })

  # Present interactive mz value plot
  output$mz_plot = ggiraph::renderGirafe({
      withProgress(message="Updating plot...", value=0, {
        x_ggr = ggiraph::girafe(ggobj=base_mz_plot(),
                                width_svg=12, height_svg=12,
                                ggiraph::opts_selection(type='single'))

        x_ggr = ggiraph::girafe_options(x_ggr, ggiraph::opts_zoom(max=5),
                                        ggiraph::opts_selection(type="multiple",
                                                                css="fill:cyan;stroke:black;opacity:0.7;"))

        return(x_ggr)
      })
  })

  # Generate plot to present tissue image (intensities module)
  output$img_mz_plot = renderImage({
    if(!is.null(input$img_user_file)){
      list(src=input$img_user_file$datapath)
    }
  }, deleteFile=F)

  # mzval PDF download
  output$download_mz_pdf = downloadHandler(
    filename = function(){
      paste0("spatial_intensity_", input$mz_selected, '_', Sys.Date(), ".pdf")
    },
    content = function(file){
      ggsave(file, plot=base_mz_plot(), device="pdf", width=8, height=6)
    }
  )

  # mzval PNG download
  output$download_mz_png = downloadHandler(
    filename = function(){
      paste0("spatial_intensity_", input$mz_selected, '_', Sys.Date(), ".png")
    },
    content = function(file){
      ggsave(file, plot=base_mz_plot(), device="png", width=8, height=6, dpi=300)
    }
  )

  # Download annotations from intensities page
  output$export_df_mz = downloadHandler(
    filename = function(){
      paste0("manual_annotations_", Sys.Date(), ".csv")
    },
    content = function(file){
      df_tmp = pixel_clusters()
      df_tmp = data.frame(id=df_tmp[['pixel_id']],
                          x_coord=df_tmp[['x_coord']],
                          y_coord=df_tmp[['y_coord']],
                          manual_from_hcdist=df_tmp[['hc_manual']],
                          manual_from_intx=df_tmp[['mz_manual']])

      write.csv(df_tmp, file, quote=T, row.names=F)
    }
  )


  ########################### DIFFERENTIAL ABUNDANCE ###########################
  diff_res = eventReactive(input$test_run_mz, {
    # Threshold of intensity to filter out metabolites
    itx_thr = as.numeric(input$user_itx_thr)

    # Subset annotations to user-selected
    df_tmp = data.frame(pixel_id=manualout_res()[['pixel_id']],
                        annots=manualout_res()[[input$annotation_test]])
    df_tmp = df_tmp[df_tmp[[2]] %in% c(input$group1_selected, input$group2_selected), ]
    df_tmp[[2]] = as.character(df_tmp[[2]])

    # Subset intensity matrix to relevant clusters
    itx_tmp = itx_res()[, df_tmp[['pixel_id']]]

    # Find metabolites with expression above user threshold
    metabs_means = sort(Matrix::rowMeans(itx_tmp), decreasing=T)
    metabs_test = names(metabs_means[ 1:ceiling(nrow(itx_tmp)*as.numeric(input$user_itx_thr)) ])

    # Subset matrix to molecules above threshold
    itx_tmp = itx_tmp[metabs_test, ]

    withProgress(message="Running tests...", value=0, {
      # Create dataframe to store DA results
      da_tests_df = data.frame(metabolite=NA, cluster=NA, avg_log_itx1=NA, avg_log_itx2=NA, log_fc=NA, test_statistic=NA, p_value=NA, adj_p_value=NA)
      if(length(unique(df_tmp[['annots']])) >= 2){
        da_tests_df = da_tests(x=itx_tmp,
                               group1=df_tmp[['pixel_id']][ df_tmp[[2]] == input$group1_selected ],
                               group2=df_tmp[['pixel_id']][ df_tmp[[2]] == input$group2_selected ],
                               group1_lbl=input$group1_selected,
                               group2_lbl=input$group2_selected,
                               method=input$diff_test, pval_adj='fdr')
      }
    })

    return(da_tests_df)
  })

  # Present DA output table
  output$da_output = DT::renderDT(
    DT::datatable(diff_res() %>%
                    dplyr::mutate(dplyr::across(c("avg_log_itx1", "avg_log_itx2", "log_fc", "p_value", "adj_p_value"), round, digits=5)) %>%
                    dplyr::rename(`Molecule`=metabolite,
                                  `Domain`=cluster,
                                  `Avg. log(intx) 1`=avg_log_itx1,
                                  `Avg. log(intx) 2`=avg_log_itx2,
                                  `log(fold change)`=log_fc,
                                  `Test statistic`=6,
                                  `p-value`=p_value,
                                  `FDR`=adj_p_value),
                  options=list(scrollX=TRUE, scrollCollapse=TRUE, pageLength=10))
  )


  ####################### SPATIAL STATISTICS CALCULATION #######################
  spatial_stats_res = eventReactive(input$test_run_spatial_stats, {
    # Number of features to calculate statistics (user-selected)
    top_features = as.integer(input$user_top_var)
    top_features = top_var()[1:top_features]

    # Extract cooordinates
    df_tmp = data.frame(pixel_id=manualout_res()[['pixel_id']],
                        x_coord=manualout_res()[['x_coord']],
                        y_coord=manualout_res()[['y_coord']])

    # Subset intensity matrix to top variable features
    itx_tmp = itx_res()[top_features, ]

    withProgress(message="Calculating spatial weights...", value=0, {
      # If too many pixels, distance matrix eats all memory. So, use KNN method instead
      k=NULL
      if(nrow(df_tmp) >= 45000){
        k = 8
      }
      spw = calculate_spatial_weights(coords=df_tmp, k=k)

    })

    withProgress(message="Calculating statistics...", value=0, {
      sp_res = spatial_stats(x=itx_tmp, spw=spw)
    })

    return(sp_res)
  })

  # Present DA output table
  output$sp_output = DT::renderDT(
    DT::datatable(spatial_stats_res() %>%
                    dplyr::mutate(dplyr::across(c("moran_i", "geary_c"), round, digits=5)) %>%
                    dplyr::rename(`Molecule`=feature,
                                  `Moran's I`=moran_i,
                                  `Geary's C`=geary_c),
                  options=list(scrollX=TRUE, scrollCollapse=TRUE, pageLength=10))
  )


  ####################### SPATIAL GRADIENTS TEST ######################
  spatial_gradients_res = eventReactive(input$test_run_spatial_gradients, {
    # Number of features to calculate statistics (user-selected)
    top_features = as.integer(input$user_top_var_grad)
    top_features = top_var()[1:top_features]

    # Extract coordinates and annotations
    df_tmp = data.frame(pixel_id=manualout_res()[['pixel_id']],
                        x_coord=manualout_res()[['x_coord']],
                        y_coord=manualout_res()[['y_coord']],
                        annots=manualout_res()[[input$annotation_test_gradients]])
    df_tmp[[4]] = as.character(df_tmp[[4]])

    # Subset intensity matrix to top variable features
    itx_tmp = itx_res()[top_features, ]

    # If too many pixels, distance matrix eats all memory. So, find neighborhoods to collapse data
    nbs_ls = NULL
    min_nb = 9
    if(nrow(df_tmp) >= 45000){
      withProgress(message="Identifying pixel neighborhoods...", value=0, {
        nbs_ls = find_adj_neighbors(coords_df=df_tmp[, 1:3])
      })
      min_nb = 1
    }

    # Perform spatial gradients tests
    spg_tests_df = data.frame(metabolite=NA, cluster=NA, avg_log_itx1=NA, avg_log_itx2=NA,
                              log_fc=NA, test_statistic=NA, p_value=NA, adj_p_value=NA)
    spg_tests_df = spatial_gradient(x=itx_tmp,
                                    sp_df=df_tmp[, 1:3],
                                    ref=df_tmp[['pixel_id']][ df_tmp[[4]] == input$ref_group ],
                                    distsumm=input$summ_type,
                                    nbs_ls=nbs_ls,
                                    min_nb=min_nb,
                                    log_dist=T)

    return(spg_tests_df)
  })

  # Present DA output table
  output$spgradient_output = DT::renderDT(
    DT::datatable(spatial_gradients_res() %>%
                    dplyr::mutate(dplyr::across(c("lm_coef", "lm_pval", "spearman_r", "spearman_r_pval", "spearman_r_pval_adj"), round, digits=5)) %>%
                    dplyr::rename(`Molecule`=molecule,
                                  `LM Slope`=lm_coef,
                                  `Slope p-value`=lm_pval,
                                  `Spearman R`=spearman_r,
                                  `Spearman p-value`=spearman_r_pval,
                                  `FDR (Spearman)`=spearman_r_pval_adj,
                                  `Comments`=pval_comment),
                  options=list(scrollX=TRUE, scrollCollapse=TRUE, pageLength=10))
  )


  ################################# BUTTONS ####################################
  # Declare reactive value to store manual ROI annotations
  # manual_annots = reactiveVal(NULL)
  # observe({
  #   manual_annots = reactiveVal(hcdist_res())
  # })

  # Button to confirm manual annotation on hcdist image
  observeEvent(input$label_confirm_hc, {
    updated_ann = pixel_clusters()
    if(input$label_input_numerichc != ""){
      cat('LABEL SELECTION - DOMAINS\n')
      user_hc_input = gsub(" ", '', input$label_input_numerichc)
      user_hc_input = unique(unlist(str_split(user_hc_input, ',')))
      rows_selected = updated_ann[['row_id']][updated_ann[['hc_manual']] %in% user_hc_input]
    } else{
      cat('LASSO SELECTION - DOMAINS\n')
      rows_selected = updated_ann[['row_id']][as.integer(input$hc_plot_selected)]
    }
    updated_ann[['hc_manual']][updated_ann[['row_id']] %in% rows_selected] = input$label_input_hc

    pixel_clusters(updated_ann)
  })

  # Button to confirm manual annotation on m/z image
  observeEvent(input$label_confirm_mz, {
    updated_ann = pixel_clusters()
    cat('LASSO SELECTION - INTENSITIES\n')
    rows_selected = updated_ann[['row_id']][as.integer(input$mz_plot_selected)]
    updated_ann[['mz_manual']][updated_ann[['row_id']] %in% rows_selected] = input$label_input_mz

    pixel_clusters(updated_ann)
  })

  # Button to run spatial statistics
  observeEvent(input$test_run_spatial_stats, {
    print("Button clicked")
  })

  # Observer to show notification if cluster file is invalid
  observeEvent(input$hcdist_user_file, {
    if (!valid_hcdist_file()){
      showNotification("Clustering file is not tab-delimited or does not contain enough columns.",
                       type="error",
                       duration=20,  # seconds
                       closeButton=TRUE
      )
    }
  })

  # Observer to show notification if intensity matrix file is invalid
  observeEvent(input$itx_user_file, {
    if (!valid_itx_file()){
      showNotification("Intensity file is not tab-delimited or does not contain enough columns.",
                       type="error",
                       duration=20,  # seconds
                       closeButton=TRUE
      )
    }
  })

  # Validation of clusters to re-annotate
  numhc_thr_v = InputValidator$new()
  numhc_thr_v$add_rule("label_input_numerichc", sv_regex("^[ ,0-9]*$", "Only integers and comma allowed"))
  numhc_thr_v$enable()

  # Validation of DA threshold value
  da_thr_v = InputValidator$new()
  da_thr_v$add_rule("user_itx_thr", sv_between(0, 1))
  da_thr_v$enable()

  # Validation of clusters to test (cant be the same cluster)
  observe({
    cl1_thr_v = InputValidator$new()
    cl1_thr_v$add_rule("group1_selected", sv_not_equal(input$group2_selected))
    cl1_thr_v$enable()

    cl2_thr_v = InputValidator$new()
    cl2_thr_v$add_rule("group2_selected", sv_not_equal(input$group1_selected))
    cl2_thr_v$enable()
  })

  # For validation of data input and enabling of modules
  input_validate = reactive({
    val = !(is.null(input$hcdist_user_file) | is.null(input$itx_user_file) )
    return(val)
  })
  output$fileUploaded = input_validate
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

  # Validate tab-delimited file for clusters
  valid_hcdist_file = reactive({
    infile = input$hcdist_user_file
    req(infile)

    first_line = readLines(infile$datapath, n=1)
    tsv_check = gregexpr('\\t', first_line)
    num_tabs = lengths(regmatches(first_line, tsv_check))

    if(num_tabs < 3){
      return(FALSE)
    }

    return(TRUE)
  })

  # Validate tab-delimited file for metabolite intensity
  valid_itx_file = reactive({
    infile = input$itx_user_file
    req(infile)

    first_line = readLines(infile$datapath, n=1)
    tsv_check = gregexpr('\\t', first_line)
    num_tabs = lengths(regmatches(first_line, tsv_check))

    if(num_tabs < 3){
      return(FALSE)
    }

    return(TRUE)
  })

  # For validation of image input
  input_imgvalidate = reactive({
    val = !is.null(input$img_user_file)
    return(val)
  })
  output$imgUploaded = input_imgvalidate
  outputOptions(output, 'imgUploaded', suspendWhenHidden=FALSE)

  # Front page documentation
  output$getting_started = renderUI({
    withMathJax({
      k = knitr::knit(input="vignettes/getting_started.Rmd", quiet=TRUE)
      HTML(markdown::markdownToHTML(k, fragment.only=T))
    })
  })

  # Spatial domains documentation
  output$spatial_domains_doc = renderUI({
    withMathJax({
      k = knitr::knit(input="vignettes/spatial_domains.Rmd", quiet=TRUE)
      HTML(markdown::markdownToHTML(k, fragment.only=T))
    })
  })

}) # CLOSE SERVER

