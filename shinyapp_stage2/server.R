library(shiny)
library(data.table)
library(tidyverse)
library(limma)
library(splines)
library(ggplot2)
library(ggrepel)
library(DelayedArray)
library(Matrix)
library(uwot)
library(spatstat.geom)
library(spatstat.explore)
library(mxfda)      
library(parallel) 
library(future)
library(promises)
library(WGCNA)
library(parallelly)
library(ggnewscale)
library(spdep)
library(sfsmisc)


plan(multisession, workers = availableCores() - 1)
options(future.globals.maxSize = 8000 * 1024^2)

# Source existing modules
source("R/de_analysis.R", local = TRUE)
source("R/volcano_plot.R", local = TRUE)
source("R/violin_plot.R", local = TRUE)

# NEW: source helper files
source("R/helpers.R", local = TRUE)
source("R/de_helpers.R", local = TRUE)
# In server.R, after loading the new helper
source("R/preprocess.R", local = TRUE)  

shinyServer(function(input, output, session) {
    
    # --------------------------------------------------------------------------
    # Data loading (unchanged)
    # --------------------------------------------------------------------------
    valid_hcdist_file <- reactive({
        infile <- input$hcdist_user_file
        req(infile)
        first_line <- readLines(infile$datapath, n = 1)
        num_tabs <- lengths(regmatches(first_line, gregexpr('\\t', first_line)))
        num_tabs >= 3
    })
    
    valid_itx_file <- reactive({
        infile <- input$itx_user_file
        req(infile)
        first_line <- readLines(infile$datapath, n = 1)
        num_tabs <- lengths(regmatches(first_line, gregexpr('\\t', first_line)))
        num_tabs >= 3
    })
    
    hcdist_res <- eventReactive(input$hcdist_user_file, {
        req(valid_hcdist_file())
        infile <- input$hcdist_user_file
        cat('READING CLUSTER TABLE\n')
        pixel_clusters <- data.table::fread(infile$datapath,
                                            check.names = FALSE, data.table = FALSE,
                                            colClasses = list(character = c("cluster")))
        hcdist_headers <- c("Pixel", "X", "Y", "tree order", "cluster", "cluster_size", "color")
        # Use clean_color from helpers
        if (sum(colnames(pixel_clusters) %in% hcdist_headers) != 7) {
            pixel_clusters <- pixel_clusters[, c(1:4)]
            colnames(pixel_clusters) <- c("pixel_id", "x_coord", "y_coord", "hc_orig")
            set.seed(12345)
            col_pal <- sample(khroma::color("smoothrainbow", force = TRUE)(length(unique(pixel_clusters[['hc_orig']]))))
            col_pal <- data.frame(hc_orig = unique(pixel_clusters[['hc_orig']]), color = col_pal)
            pixel_clusters <- pixel_clusters %>% dplyr::left_join(., col_pal, by = 'hc_orig')
            pixel_clusters$color <- sapply(pixel_clusters$color, clean_color)
        } else if (sum(colnames(pixel_clusters) %in% hcdist_headers) == 7) {
            pixel_clusters <- pixel_clusters %>%
                dplyr::select("pixel_id" = "Pixel", 'x_coord' = "X", 'y_coord' = "Y",
                              'hc_orig' = "cluster", "color")
            pixel_clusters$color <- sapply(pixel_clusters$color, clean_color)
        }
        pixel_clusters[['hc_orig']] <- as.character(pixel_clusters[['hc_orig']])
        pixel_clusters[['x_coord']] <- as.integer(pixel_clusters[['x_coord']])
        pixel_clusters[['y_coord']] <- as.integer(pixel_clusters[['y_coord']])
        pixel_clusters[['hc_manual']] <- pixel_clusters[['hc_orig']]
        pixel_clusters[['mz_manual']] <- rep('Default', nrow(pixel_clusters))
        pixel_clusters <- pixel_clusters %>%
            tibble::add_column(row_id = as.integer(c(1:nrow(.))), .before = 1)
        return(pixel_clusters)
    })
    
    itx_in <- eventReactive(input$itx_user_file, {
        req(valid_itx_file())
        withProgress(message = "Reading intensity data...", value = 0, {
            infile <- input$itx_user_file
            intx <- vroom::vroom(infile$datapath, guess_max = 10, show_col_types = FALSE)
            intx <- tibble::column_to_rownames(intx, var = colnames(intx)[1])
            if (prod(dim(intx)) < 1e6) {
                intx <- as.matrix(intx)
            } else {
                intx <- DelayedArray(as.matrix(intx))
            }
        })
        return(intx)
    })
    
    pixel_clusters <- reactiveVal()
    observe({ pixel_clusters(hcdist_res()) })
    
    itx_res <- reactiveVal()
    observe({ itx_res(itx_in()) })
    
    neighbors <- reactiveVal(NULL)
    
    # --------------------------------------------------------------------------
    # Manual annotation (ROI) – unchanged
    # --------------------------------------------------------------------------
    hc_plot_data <- reactive({
        dat_tmp <- pixel_clusters()
        req(dat_tmp)
        filter_val <- input$cluster_filter
        if (!is.null(filter_val) && filter_val != "All") {
            dat_tmp <- dat_tmp %>% dplyr::filter(hc_manual == filter_val)
        }
        dat_tmp
    })
    
    selected_hc_pixels <- reactive({
        req(input$hc_plot_selected, hc_plot_data())
        idx <- as.integer(input$hc_plot_selected)
        idx <- idx[idx > 0 & idx <= nrow(hc_plot_data())]
        hc_plot_data()$pixel_id[idx]
    })
    
    mz_plot_data_for_selection <- reactive({
        req(input$mz_selected, itx_in(), pixel_clusters())
        mz_val <- input$mz_selected
        dat_base <- as.data.frame(pixel_clusters())
        int_vec <- itx_in()[mz_val, , drop = TRUE]
        plot_df <- data.frame(
            pixel_id = names(int_vec),
            mz_intx  = as.numeric(int_vec),
            stringsAsFactors = FALSE
        )
        plot_df <- merge(dat_base, plot_df, by = "pixel_id", all.x = FALSE, all.y = FALSE)
        plot_df <- plot_df[, c("pixel_id", "x_coord", "y_coord", "mz_intx")]
        colnames(plot_df) <- c("Pixel", "X", "Y", "mz_intx")
        plot_df$mz_intx[is.na(plot_df$mz_intx)] <- 0
        plot_df
    })
    
    selected_mz_pixels <- reactive({
        req(input$mz_plot_selected, mz_plot_data_for_selection())
        idx <- as.integer(input$mz_plot_selected)
        idx <- idx[idx > 0 & idx <= nrow(mz_plot_data_for_selection())]
        mz_plot_data_for_selection()$Pixel[idx]
    })
    
    observeEvent(input$label_confirm_hc, {
        req(input_validate())
        cl <- pixel_clusters()
        new_label <- input$label_input_hc
        if (is.null(new_label) || new_label == "") {
            showNotification("Please enter a label name.", type = "warning")
            return()
        }
        if (!is.null(input$label_input_numerichc) && input$label_input_numerichc != "") {
            cat("Re‑annotating clusters from numeric input\n")
            user_hc_input <- gsub(" ", "", input$label_input_numerichc)
            clusters_to_change <- unique(unlist(strsplit(user_hc_input, ",")))
            cl$hc_manual[cl$hc_manual %in% clusters_to_change] <- new_label
        } else {
            pixel_ids <- selected_hc_pixels()
            if (length(pixel_ids) == 0) {
                showNotification("No pixels selected. Use the lasso tool.", type = "warning")
                return()
            }
            cl$hc_manual[cl$pixel_id %in% pixel_ids] <- new_label
        }
        pixel_clusters(cl)
        showNotification(paste("Annotation updated to", new_label), type = "message")
    })
    
    observeEvent(input$label_confirm_mz, {
        req(input_validate())
        cl <- pixel_clusters()
        new_label <- input$label_input_mz
        if (is.null(new_label) || new_label == "") {
            showNotification("Please enter a label name.", type = "warning")
            return()
        }
        pixel_ids <- selected_mz_pixels()
        if (length(pixel_ids) == 0) {
            showNotification("No pixels selected. Use the lasso tool.", type = "warning")
            return()
        }
        cl$mz_manual[cl$pixel_id %in% pixel_ids] <- new_label
        pixel_clusters(cl)
        showNotification(paste("Annotation updated to", new_label), type = "message")
    })
    
    # --------------------------------------------------------------------------
    # Input validation & UI helpers (unchanged)
    # --------------------------------------------------------------------------
    input_validate <- reactive({
        !(is.null(input$hcdist_user_file) | is.null(input$itx_user_file))
    })
    output$fileUploaded <- input_validate
    outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
    
    output$getting_started <- renderUI({
        includeMarkdown("getting_started.md")
    })
    
    observe({
        req(pixel_clusters())
        clusters <- sort(unique(pixel_clusters()$hc_manual))
        updateSelectInput(session, "cluster_filter",
                          choices = c("All" = "All", as.character(clusters)),
                          selected = input$cluster_filter %||% "All")
    })
    observeEvent(input$label_confirm_hc, {
        req(pixel_clusters())
        clusters <- sort(unique(pixel_clusters()$hc_manual))
        updateSelectInput(session, "cluster_filter",
                          choices = c("All" = "All", as.character(clusters)),
                          selected = input$cluster_filter %||% "All")
    }, ignoreInit = TRUE)
    
    # --------------------------------------------------------------------------
    # Pixel collapsing (unchanged)
    # --------------------------------------------------------------------------
    observeEvent(list(input$hcdist_user_file, input$itx_user_file, input$collapse_par), {
        req(input$hcdist_user_file, input$itx_user_file, valid_hcdist_file())
        min_nb <- as.integer(input$collapse_par)
        hc_out <- pixel_clusters()
        coords_df <- hc_out[, c('pixel_id', 'x_coord', 'y_coord')]
        systematic_sub <- FALSE
        if (nrow(coords_df) >= 50000) systematic_sub <- TRUE
        if (nrow(coords_df) >= 20e6 | !is.integer(min_nb) | length(min_nb) == 0) {
            min_nb <- 21
        } else if (systematic_sub) {
            min_nb <- NA_integer_
        }
        if (!is.na(min_nb)) {
            cat('COLLAPSING PIXELS (parallel)\n')
            withProgress(message = "Collapsing pixels", value = 0, {
                nbs_ls <- find_adj_neighbors(coords_df = coords_df, n = min_nb)
                incProgress(0.3, detail = "Summarizing abundances in parallel")
                itx_mtx <- itx_in()
                library(future.apply)
                collapse_itx <- future_lapply(nbs_ls, function(i) {
                    itx_tmp <- DelayedMatrixStats::rowMedians(itx_mtx[, i, drop = FALSE])
                    as.matrix(itx_tmp)
                }, future.seed = TRUE) %>% do.call('cbind', .)
                colnames(collapse_itx) <- names(nbs_ls)
                incProgress(0.6, detail = "Summarizing labels in parallel")
                collapse_labs <- future_lapply(nbs_ls, function(i) {
                    lab_tmp <- get_most_common_category(i, hc_out)
                    as.data.frame(lab_tmp)
                }, future.seed = TRUE) %>%
                    do.call('rbind', .) %>%
                    rownames_to_column('pixel_id') %>%
                    dplyr::rename(pixel_id = 1, hc_orig = 2)
                incProgress(1)
            })
            neighbors(list(nbs_ls, collapse_itx, collapse_labs))
        } else if (systematic_sub) {
            cat('SYSTEMATIC SUBSAMPLING\n')
            subsampled_rows <- hc_out[seq(1, nrow(hc_out), by = 10), ]
            pixel_clusters(subsampled_rows)
            itx_mtx <- itx_in()
            pix_to_keep <- intersect(subsampled_rows$pixel_id, colnames(itx_mtx))
            if (length(pix_to_keep) == 0) {
                showNotification("No matching pixel IDs", type = "error")
                return()
            }
            itx_res(itx_mtx[, pix_to_keep, drop = FALSE])
            neighbors(list(NULL, NULL, NULL))
        } else {
            neighbors(list(NULL, NULL, NULL))
        }
    }, ignoreInit = TRUE)
    
    # --------------------------------------------------------------------------
    # Most variable features (unchanged)
    # --------------------------------------------------------------------------
    top_var <- reactive({
        cat('IDENTIFYING MOST VARIABLE FEATURES (fast)\n')
        nbs_ls <- neighbors()
        if (is.null(nbs_ls[[1]])) {
            itx_tmp <- itx_res()
        } else {
            itx_tmp <- nbs_ls[[2]]
        }
        if (!is.null(itx_tmp)) {
            if (inherits(itx_tmp, "DelayedArray")) {
                sds <- DelayedMatrixStats::rowSds(itx_tmp)
            } else {
                sds <- matrixStats::rowSds(itx_tmp)
            }
            var_mz <- sort(sds, decreasing = TRUE)
            names(var_mz) <- rownames(itx_tmp)
            return(names(var_mz))
        }
    })
    
    observe({
        updateSelectizeInput(session, 'mz_selected', choices = top_var(), server = TRUE)
    })
    
    # ---- Differential Abundance group selectors (unchanged) ----
    observe({
        req(input_validate(), pixel_clusters(), input$annotation_test)
        cluster_opts <- unique(pixel_clusters()[[input$annotation_test]])
        cluster_opts <- cluster_opts[!is.na(cluster_opts)]
        updateSelectizeInput(session, 'group1_selected',
                             choices = cluster_opts,
                             server = TRUE,
                             selected = if (length(cluster_opts) > 0) cluster_opts[1] else NULL)
    })
    observe({
        req(input_validate(), pixel_clusters(), input$annotation_test)
        cluster_opts <- unique(pixel_clusters()[[input$annotation_test]])
        cluster_opts <- cluster_opts[!is.na(cluster_opts)]
        updateSelectizeInput(session, 'group2_selected',
                             choices = cluster_opts,
                             server = TRUE,
                             selected = if (length(cluster_opts) > 1) cluster_opts[2] else cluster_opts[1])
    })
    
    # ---- Spatial gradients reference cluster (unchanged) ----
    observe({
        req(input_validate(), pixel_clusters(), input$annotation_test_gradients)
        cluster_opts <- unique(pixel_clusters()[[input$annotation_test_gradients]])
        cluster_opts <- cluster_opts[!is.na(cluster_opts)]
        updateSelectizeInput(session, 'ref_group',
                             choices = cluster_opts,
                             server = TRUE,
                             selected = if (length(cluster_opts) >= 3) cluster_opts[3] else cluster_opts[1])
    })
    
    
    # =============================================================================
    # PREPROCESSING MODULE
    # =============================================================================
    
    # Source the preprocessing orchestrator
    source("R/preprocess.R", local = TRUE)
    
    # -----------------------------------------------------------------------------
    # Reactive values
    # -----------------------------------------------------------------------------
    dataLoaded <- reactiveVal(FALSE)          # controls sidebar visibility
    preproc_status <- reactiveVal("Waiting for input")
    preproc_results <- reactiveVal(NULL)      # stores output paths (for downloads)
    generated_data <- reactiveVal(NULL)       # stores parsed data after preprocessing
    preproc_done <- reactiveVal(FALSE)        # flag to show load button
    
    # Store the path to the extracted .imzML file after upload
    if (is.null(session$userData$imzml_path)) {
        session$userData$imzml_path <- reactiveVal(NULL)
    }
    
    # -----------------------------------------------------------------------------
    # File upload handler: accept imzML+ibd or ZIP
    # -----------------------------------------------------------------------------
    observeEvent(input$raw_imzml_files, {
        req(input$raw_imzml_files)
        
        files <- input$raw_imzml_files
        temp_dir <- tempdir()
        
        # Check if it's a ZIP (single file with .zip extension)
        if (length(files) == 1 && grepl("\\.zip$", files$name[1], ignore.case = TRUE)) {
            # Extract ZIP
            zip_path <- files$datapath[1]
            tryCatch({
                unzip(zip_path, exdir = temp_dir)
            }, error = function(e) {
                showNotification("Failed to extract ZIP. Please check the file.", type = "error")
                return()
            })
            imzml_files <- list.files(temp_dir, pattern = "\\.imzML$", ignore.case = TRUE, full.names = TRUE)
            if (length(imzml_files) == 0) {
                showNotification("No .imzML file found in ZIP.", type = "error")
                return()
            }
            session$userData$imzml_path(imzml_files[1])
            preproc_status(paste("ZIP extracted. Found:", basename(imzml_files[1])))
            showNotification("ZIP ready for preprocessing.", type = "message")
            
        } else {
            # Direct upload of imzML and ibd (multiple files)
            # Find the .imzML file among the uploaded files
            imzml_idx <- grep("\\.imzML$", files$name, ignore.case = TRUE)
            if (length(imzml_idx) == 0) {
                showNotification("No .imzML file found among uploaded files.", type = "error")
                return()
            }
            imzml_file <- files$datapath[imzml_idx[1]]
            imzml_orig_name <- files$name[imzml_idx[1]]
            
            # Find corresponding .ibd (or .IBD)
            ibd_idx <- grep("\\.ibd$", files$name, ignore.case = TRUE)
            if (length(ibd_idx) == 0) {
                showNotification("No .ibd file found. The .imzML requires its .ibd.", type = "error")
                return()
            }
            ibd_file <- files$datapath[ibd_idx[1]]
            ibd_orig_name <- files$name[ibd_idx[1]]
            
            # Copy both files to a common temporary directory with their original names
            target_dir <- file.path(temp_dir, "imzml_upload")
            dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)
            
            file.copy(imzml_file, file.path(target_dir, imzml_orig_name), overwrite = TRUE)
            file.copy(ibd_file, file.path(target_dir, ibd_orig_name), overwrite = TRUE)
            
            imzml_path <- file.path(target_dir, imzml_orig_name)
            session$userData$imzml_path(imzml_path)
            preproc_status(paste("Uploaded:", basename(imzml_path)))
            showNotification("Files ready for preprocessing.", type = "message")
        }
    })
    
    # -----------------------------------------------------------------------------
    # Main preprocessing observer
    # -----------------------------------------------------------------------------
    observeEvent(input$run_preprocess, {
        # Ensure we have a valid imzML path
        imzml_path <- session$userData$imzml_path()
        if (is.null(imzml_path)) {
            showNotification("Please upload .imzML and .ibd files (or a ZIP archive) first.", type = "warning")
            return()
        }
        
        # Collect all parameters from UI (including advanced options)
        params <- list(
            cores = input$prep_cores,
            n_mz_clusters = input$prep_n_mz_clusters,
            filter_present = input$prep_filter_present,
            filter_unlog_sd = input$prep_filter_unlog_sd,
            filter_log2_sd = input$prep_filter_log2_sd,
            minkowski_p = input$prep_minkowski_p,
            distpow = input$prep_distpow,
            linkage = input$prep_linkage,
            log2 = input$prep_log2,
            mean_center = input$prep_mean_center,
            unit_variance = input$prep_unit_variance,
            transpose = input$prep_transpose
        )
        
        # Create a unique output directory
        out_dir <- file.path(tempdir(), paste0("hcdist_run_", Sys.getpid()))
        
        # Reset flags and update status
        preproc_done(FALSE)
        generated_data(NULL)
        dataLoaded(FALSE)
        preproc_status("Running hcdist pipeline... (this may take a while)")
        shinyjs::disable("run_preprocess")
        
        # Launch the pipeline in a background future
        future({
            preprocess_raw_data(imzml_path, out_dir, params, hcdist_root = "../hcdist_stage1")
        }) %...>% (function(res) {
            # ---- Success: read output files ----
            
            # Read cluster file
            cluster_df <- read.table(res$cluster_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            required <- c("Pixel", "X", "Y", "cluster")
            if (!all(required %in% colnames(cluster_df))) {
                if (ncol(cluster_df) >= 4) {
                    colnames(cluster_df)[1:4] <- required
                } else {
                    stop("Cluster file format incorrect.")
                }
            }
            # Build internal data frame
            cluster_df <- cluster_df[, required]
            colnames(cluster_df) <- c("pixel_id", "x_coord", "y_coord", "hc_orig")
            cluster_df$hc_orig <- as.character(cluster_df$hc_orig)   # <-- ensure character
            cluster_df$hc_manual <- cluster_df$hc_orig               # <-- already character
            cluster_df$mz_manual <- "Default"
            cluster_df$row_id <- seq_len(nrow(cluster_df))
            # Assign colors based on unique clusters
            uniq_clusters <- unique(cluster_df$hc_orig)
            color_pal <- rainbow(length(uniq_clusters))
            names(color_pal) <- uniq_clusters
            cluster_df$color <- color_pal[cluster_df$hc_orig]        # <-- now works with character
            # Reorder to match app's expected columns
            cluster_df <- cluster_df[, c("row_id", "pixel_id", "x_coord", "y_coord", 
                                         "hc_orig", "hc_manual", "mz_manual", "color")]
            
            # Read intensity matrix
            int_mat <- read.table(res$intensity_file, header = TRUE, sep = "\t", row.names = 1)
            int_mat <- as.matrix(int_mat)
            
            # ---- Store generated data (not loaded yet) ----
            generated_data(list(
                cluster_df = cluster_df,
                int_mat = int_mat,
                cluster_file = res$cluster_file,
                intensity_file = res$intensity_file
            ))
            
            # Update status and enable button
            preproc_status("Preprocessing completed! Click 'Load generated data' to use it.")
            shinyjs::enable("run_preprocess")
            
            # Show the "Load" button and download buttons
            preproc_done(TRUE)
            
            # Store results for downloads
            preproc_results(res)
            
            showNotification("Preprocessing finished. Use the 'Load generated data' button to import.", type = "message")
            
        }) %...!% (function(err) {
            preproc_status(paste("Error:", err$message))
            shinyjs::enable("run_preprocess")
            showNotification(paste("Preprocessing failed:", err$message), type = "error")
        })
    })
    
    # -----------------------------------------------------------------------------
    # Observer for "Load generated data" button
    # -----------------------------------------------------------------------------
    observeEvent(input$load_generated_data, {
        req(generated_data())
        data <- generated_data()
        
        # Assign to the app's main reactive values
        pixel_clusters(data$cluster_df)
        itx_res(data$int_mat)
        
        # Set dataLoaded to TRUE so sidebar appears
        dataLoaded(TRUE)
        
        # Hide the load button (optional)
        preproc_done(FALSE)
        preproc_status("Data loaded successfully!")
        
        showNotification("Generated data loaded into app.", type = "message")
    })
    
    # -----------------------------------------------------------------------------
    # Outputs for the preprocessing UI
    # -----------------------------------------------------------------------------
    
    # Status text
    output$preprocess_status <- renderPrint({
        cat(preproc_status())
    })
    
    # Preview of generated cluster data (first 10 rows)
    output$prep_preview <- renderTable({
        req(generated_data())
        head(generated_data()$cluster_df, 10)
    })
    
    # Conditional: show download buttons and load button after preprocessing
    output$preprocessing_done <- reactive({
        preproc_done()
    })
    outputOptions(output, "preprocessing_done", suspendWhenHidden = FALSE)
    
    # Conditional: show analysis tabs after data is loaded
    output$dataLoaded <- reactive({
        dataLoaded()
    })
    outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
    
    # -----------------------------------------------------------------------------
    # Download handlers for the generated files
    # -----------------------------------------------------------------------------
    
    output$download_prep_cluster <- downloadHandler(
        filename = function() paste0("hcdist_result_", Sys.Date(), ".txt"),
        content = function(file) {
            req(preproc_results())
            file.copy(preproc_results()$cluster_file, file)
        }
    )
    
    output$download_prep_intensity <- downloadHandler(
        filename = function() paste0("intensity_table_", Sys.Date(), ".tsv"),
        content = function(file) {
            req(preproc_results())
            file.copy(preproc_results()$intensity_file, file)
        }
    )
    # --------------------------------------------------------------------------
    # SPATIAL DOMAINS PLOT 
    # --------------------------------------------------------------------------
    base_hc_plot <- reactive({
        dat_tmp <- pixel_clusters()
        req(dat_tmp)
        
        filter_val <- input$cluster_filter
        if (!is.null(filter_val) && filter_val != "All") {
            dat_tmp <- dat_tmp %>% dplyr::filter(hc_manual == filter_val)
        }
        
        if (nrow(dat_tmp) == 0) {
            return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                       label = paste("No pixels for cluster", filter_val)) + theme_void())
        }
        
        # Convert to factor to ensure discrete scale (fix "Continuous value supplied to discrete scale")
        labs <- dat_tmp[['hc_manual']]
        alpha_val <- as.numeric(input$alpha_value_hc)
        p <- make_plot_hcdist(pixel_cl = dat_tmp, labels = labs, spot_alpha = alpha_val)
        
        if (!is.null(input$white_bg_hc) && input$white_bg_hc) {
            p <- p + theme(panel.background = element_rect(fill = "white", colour = NA),
                           plot.background = element_rect(fill = "white", colour = NA))
        }
        p
    })
    output$hc_plot <- renderPlot({ base_hc_plot() }, height = 600, width = 600, res = 96)
    
    output$img_hc_plot <- renderImage({
        if (!is.null(input$img_user_file)) list(src = input$img_user_file$datapath)
    }, deleteFile = FALSE)
    
    output$download_hc_pdf <- downloadHandler(
        filename = function() paste0("spatial_domains_", Sys.Date(), ".pdf"),
        content = function(file) ggsave(file, plot = base_hc_plot(), device = "pdf", width = 8, height = 6)
    )
    output$download_hc_png <- downloadHandler(
        filename = function() paste0("spatial_domains_", Sys.Date(), ".png"),
        content = function(file) ggsave(file, plot = base_hc_plot(), device = "png", width = 8, height = 6, dpi = 300)
    )
    output$export_df_hc <- downloadHandler(
        filename = function() paste0("manual_annotations_", Sys.Date(), ".csv"),
        content = function(file) {
            df_tmp <- pixel_clusters()
            df_tmp <- data.frame(
                id = df_tmp[['pixel_id']],
                x_coord = df_tmp[['x_coord']],
                y_coord = df_tmp[['y_coord']],
                manual_from_hcdist = df_tmp[['hc_manual']],
                manual_from_intx = df_tmp[['mz_manual']]
            )
            write.csv(df_tmp, file, quote = TRUE, row.names = FALSE)
        }
    )
    
    # --------------------------------------------------------------------------
    # INTENSITIES PLOT (unchanged, but uses helpers)
    # --------------------------------------------------------------------------
    mz_plot_settings <- eventReactive(
        input$apply_mz_settings,
        {
            list(
                alpha     = input$alpha_value_mz,
                color     = input$mz_spot_color,
                white_bg  = input$white_bg_mz,
                contrast  = input$mz_contrast,
                point_size = input$mz_point_size
            )
        },
        ignoreInit = TRUE
    )
    
    contrast_input <- shiny::debounce(reactive(input$mz_contrast), 450)
    
    
    mz_plot_data <- reactive({
        req(input$mz_selected, itx_res(), pixel_clusters())
        mz_val <- input$mz_selected
        dat_base <- as.data.frame(pixel_clusters())
        int_vec <- itx_res()[mz_val, , drop = TRUE]
        plot_df <- data.frame(
            pixel_id = names(int_vec),
            mz_intx  = as.numeric(int_vec),
            stringsAsFactors = FALSE
        )
        plot_df <- merge(dat_base, plot_df, by = "pixel_id", all.x = FALSE, all.y = FALSE)
        plot_df <- plot_df[, c("x_coord", "y_coord", "mz_intx")]
        plot_df$mz_intx[is.na(plot_df$mz_intx)] <- 0
        plot_df
    })
    
    base_mz_plot <- reactive({
        req(mz_plot_data())
        dat <- mz_plot_data()
        settings <- mz_plot_settings()
        
        contrast_val <- isolate(contrast_input())
        contrast_val <- max(0.3, min(5, suppressWarnings(as.numeric(contrast_val %||% 1))))
        gamma <- 1 / contrast_val
        dat$mz_intx <- dat$mz_intx ^ gamma
        
        rng <- range(dat$mz_intx, na.rm = TRUE)
        if (diff(rng) > 0) {
            dat$mz_intx <- (dat$mz_intx - rng[1]) / diff(rng)
        } else {
            dat$mz_intx <- 0
        }
        
        spot_color <- settings$color %||% "#1f78b4"
        point_size <- as.numeric(settings$point_size %||% 3)
        alpha_val <- as.numeric(settings$alpha %||% 0.7)
        
        p <- ggplot(dat, aes(x = x_coord, y = y_coord, color = mz_intx)) +
            geom_point(size = point_size, alpha = alpha_val) +
            scale_color_gradientn(
                colours = colorRampPalette(c("white", spot_color))(256),
                name = "Intensity",
                limits = c(0, 1)
            ) +
            scale_y_reverse() +
            coord_fixed() +
            labs(title = input$mz_selected) +
            theme_void() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 8)
            )
        
        if (isTRUE(settings$white_bg)) {
            p <- p + theme(
                panel.background = element_rect(fill = "white", colour = NA),
                plot.background = element_rect(fill = "white", colour = NA),
                legend.background = element_rect(fill = "white"),
                legend.text = element_text(colour = "black"),
                legend.title = element_text(colour = "black"),
                plot.title = element_text(colour = "black")
            )
        } else {
            p <- p + theme(
                panel.background = element_rect(fill = "black", colour = NA),
                plot.background = element_rect(fill = "black", colour = NA),
                legend.background = element_rect(fill = "black"),
                legend.text = element_text(colour = "white"),
                legend.title = element_text(colour = "white"),
                plot.title = element_text(colour = "white")
            )
        }
        p
    })
    
    output$mz_plot <- ggiraph::renderGirafe({
        req(base_mz_plot())
        ggiraph::girafe(
            ggobj = base_mz_plot(),
            width_svg = 14,
            height_svg = 14,
            pointsize = 12,
            options = list(
                ggiraph::opts_zoom(min = 0.4, max = 10),
                ggiraph::opts_selection(type = "multiple", css = "fill:cyan;stroke:black;opacity:0.9;")
            )
        )
    })
    
    output$download_mz_pdf <- downloadHandler(
        filename = function() paste0("spatial_intensity_", input$mz_selected, "_", Sys.Date(), ".pdf"),
        content = function(file) {
            ggsave(file, plot = base_mz_plot(), device = "pdf", width = 12, height = 10, dpi = 600)
        }
    )
    output$download_mz_png <- downloadHandler(
        filename = function() paste0("spatial_intensity_", input$mz_selected, "_", Sys.Date(), ".png"),
        content = function(file) {
            ggsave(file, plot = base_mz_plot(), device = "png", width = 12, height = 10, dpi = 600)
        }
    )
    
    output$img_mz_plot <- renderImage({
        if (!is.null(input$img_user_file)) {
            list(src = input$img_user_file$datapath)
        }
    }, deleteFile = FALSE)
    
    output$export_df_mz <- downloadHandler(
        filename = function() paste0("manual_annotations_", Sys.Date(), ".csv"),
        content = function(file) {
            df_tmp <- pixel_clusters()
            df_tmp <- data.frame(
                id = df_tmp$pixel_id,
                x_coord = df_tmp$x_coord,
                y_coord = df_tmp$y_coord,
                manual_from_hcdist = df_tmp$hc_manual,
                manual_from_intx = df_tmp$mz_manual
            )
            write.csv(df_tmp, file, quote = TRUE, row.names = FALSE)
        }
    )
    
    observeEvent(input$force_gc, {
        gc(full = TRUE)
        showNotification("Memory cleaned", type = "message")
    })
    
    # ----------------------------------------------------------------------------
    # Differential Abundance (now supports limma with "vs All Others")
    # ----------------------------------------------------------------------------
    
    de_results <- reactiveVal(NULL)
    
    observeEvent(input$test_run_mz, {
        req(input_validate())
        req(input$group1_selected, input$diff_test)
        if (!input$compare_to_all) req(input$group2_selected)
        
        showModal(modalDialog("Running differential analysis...", footer = NULL))
        
        tryCatch({
            withProgress(message = "Running differential analysis...", value = 0, {
                nbs_data <- neighbors()
                incProgress(0.05, detail = "Preparing data")
                
                if (!is.null(nbs_data) && !is.null(nbs_data[[2]]) && !is.null(nbs_data[[3]])) {
                    itx_use <- nbs_data[[2]]
                    cluster_df_use <- nbs_data[[3]]
                    if (!"hc_manual" %in% colnames(cluster_df_use)) {
                        cluster_df_use$hc_manual <- cluster_df_use$hc_orig
                    }
                    cat("DE using collapsed data (", ncol(itx_use), "blocks)\n")
                } else {
                    itx_use <- itx_res()
                    cluster_df_use <- pixel_clusters()
                    cat("DE using original data (", ncol(itx_use), "pixels)\n")
                }
                
                incProgress(0.1, detail = "Transforming intensities")
                itx_trans <- apply_intensity_transform(itx_use, input$intensity_transform)
                
                de_annot_col <- input$annotation_test
                de_progress <- function(step, detail = NULL) {
                    total_range <- 0.7
                    incProgress(amount = step * total_range, detail = detail)
                }
                
                incProgress(0.25, detail = "Running DE test")
                
                # ---- Handle limma separately because it needs coordinates ----
                if (input$diff_test == "limma_spatial") {
                    # Determine group2: if compare_to_all is TRUE, use "All_Others"
                    group2 <- if (input$compare_to_all) "All_Others" else input$group2_selected
                    res <- run_de_analysis(
                        itx = itx_trans,
                        cluster_df = cluster_df_use,
                        cluster_col = de_annot_col,
                        group1 = input$group1_selected,
                        group2 = group2,
                        method = "limma_spatial",
                        pval_adj = "fdr",
                        progress_callback = de_progress
                    )
                } else {
                    # Existing code for all other methods
                    if (input$compare_to_all) {
                        res <- run_de_vs_all(
                            itx = itx_trans,
                            cluster_df = cluster_df_use,
                            cluster_col = de_annot_col,
                            group1 = input$group1_selected,
                            method = input$diff_test,
                            pval_adj = "fdr",
                            progress_callback = de_progress
                        )
                    } else {
                        res <- run_de_analysis_parallel(
                            itx = itx_trans,
                            cluster_df = cluster_df_use,
                            cluster_col = de_annot_col,
                            group1 = input$group1_selected,
                            group2 = input$group2_selected,
                            method = input$diff_test,
                            pval_adj = "fdr",
                            progress_callback = de_progress
                        )
                    }
                }
                
                incProgress(0.95, detail = "Finalizing results")
                de_results(res)
                
                output$da_output <- DT::renderDT({
                    DT::datatable(res, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
                })
                
                incProgress(1, detail = "Done")
                showNotification("Differential analysis completed successfully.", type = "message")
            })
        }, error = function(e) {
            showNotification(paste("Error in DE analysis:", e$message), type = "error")
            cat("DE error:", e$message, "\n")
        }, finally = {
            removeModal()
        })
    })
    
    # ---- Volcano plot (uses helper get_cluster_colors? Not needed) ----
    volcano_gg <- reactive({
        req(de_results())
        res <- de_results()
        
        stat_col <- if ("test_statistic" %in% colnames(res)) "test_statistic" else "statistic"
        p_col <- if ("p_value" %in% colnames(res)) "p_value" else "P.Value"
        fdr_col <- if ("adj_p_value" %in% colnames(res)) "adj_p_value" else "adj.P.Val"
        fc_cutoff <- input$volcano_fc_cutoff
        
        if (all(is.na(res[[p_col]]))) {
            res$sig_group <- "Not significant"
        } else {
            res <- res %>%
                dplyr::mutate(
                    sig_group = dplyr::case_when(
                        abs(log_fc) >= fc_cutoff & .data[[fdr_col]] < input$volcano_fdr_cutoff & log_fc > 0 ~ "Up",
                        abs(log_fc) >= fc_cutoff & .data[[fdr_col]] < input$volcano_fdr_cutoff & log_fc < 0 ~ "Down",
                        TRUE ~ "Not significant"
                    )
                )
        }
        
        res <- res %>%
            dplyr::mutate(
                y_value = if (input$volcano_y == "stat") abs(.data[[stat_col]]) else -log10(.data[[p_col]] + 1e-300)
            )
        res$feature <- res$metabolite
        res <- res[!is.na(res$y_value) & is.finite(res$y_value), ]
        if (nrow(res) == 0) return(NULL)
        
        top_features <- res %>%
            dplyr::filter(sig_group != "Not significant") %>%
            dplyr::arrange(desc(y_value)) %>%
            head(input$volcano_top_n)
        
        if (!is.null(input$volcano_extra) && input$volcano_extra != "") {
            extra_names <- trimws(unlist(strsplit(input$volcano_extra, ",")))
            extra <- res %>% dplyr::filter(feature %in% extra_names)
            top_features <- dplyr::bind_rows(top_features, extra) %>% distinct()
        }
        
        y_label <- if (input$volcano_y == "stat") "Test statistic" else expression(-log[10]("p-value"))
        x_label <- if (input$intensity_transform == "log2") expression(log[2] ~ "Fold Change") else "Raw fold change (difference of means)"
        
        p <- ggplot(res, aes(x = log_fc, y = y_value)) +
            geom_point(aes(color = sig_group), alpha = 0.6, size = 1.5) +
            scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "grey70")) +
            guides(color = guide_legend(title = NULL)) +
            geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40") +
            ggrepel::geom_text_repel(data = top_features, aes(label = feature),
                                     size = input$volcano_label_size, max.overlaps = Inf) +
            theme_bw() +
            labs(x = x_label, y = y_label,
                 title = paste("Volcano plot:", input$group1_selected, "vs", 
                               if(input$compare_to_all) "All Others" else input$group2_selected)) +
            theme(axis.title = element_text(size = input$volcano_axis_title_size),
                  axis.text = element_text(size = input$volcano_axis_text_size),
                  plot.title = element_text(size = input$volcano_title_size, hjust = 0.5),
                  legend.text = element_text(size = input$volcano_legend_text_size))
        
        if (input$volcano_y == "neglogp") {
            p <- p + geom_hline(yintercept = -log10(input$volcano_fdr_cutoff), linetype = "dashed", color = "grey40")
        }
        return(p)
    })
    
    output$volcano_plot <- plotly::renderPlotly({
        req(volcano_gg())
        plotly::ggplotly(volcano_gg(), tooltip = c("x", "y", "colour", "label")) %>%
            plotly::layout(dragmode = "lasso")
    })
    
    output$download_volcano_png <- downloadHandler(
        filename = function() paste0("volcano_plot_", Sys.Date(), ".png"),
        content = function(file) {
            ggsave(file, plot = volcano_gg(), device = "png", width = 8, height = 6, dpi = 300)
        }
    )
    output$download_volcano_pdf <- downloadHandler(
        filename = function() paste0("volcano_plot_", Sys.Date(), ".pdf"),
        content = function(file) {
            ggsave(file, plot = volcano_gg(), device = "pdf", width = 8, height = 6)
        }
    )
    
    # ----------------------------------------------------------------------------
    # Violin plot (uses get_cluster_colors and apply_intensity_transform helpers)
    # ----------------------------------------------------------------------------
    
    violin_df_single <- reactive({
        req(input$violin_feature)
        df <- pixel_clusters()
        group_col <- "hc_manual"
        
        if (!is.null(input$violin_clusters) && length(input$violin_clusters) > 0) {
            df <- df[df[[group_col]] %in% input$violin_clusters, ]
        }
        req(nrow(df) > 0)
        
        vals <- as.numeric(itx_res()[input$violin_feature, df$pixel_id, drop = TRUE])
        vals <- apply_intensity_transform(vals, input$violin_transform)
        data.frame(cluster = df[[group_col]], intensity = vals, stringsAsFactors = FALSE)
    })
    
    violin_df_all <- reactive({
        df_pixel <- pixel_clusters()
        group_col <- "hc_manual"
        
        if (!is.null(input$violin_clusters) && length(input$violin_clusters) > 0) {
            df_pixel <- df_pixel[df_pixel[[group_col]] %in% input$violin_clusters, ]
        }
        req(nrow(df_pixel) > 0)
        
        clusters <- unique(df_pixel[[group_col]])
        itx_mat <- as.matrix(itx_res()[, df_pixel$pixel_id])
        if (input$violin_transform == "log2") {
            itx_mat <- log2(itx_mat + 1)
        }
        all_intensities <- list()
        for (cl in clusters) {
            pix <- df_pixel$pixel_id[df_pixel[[group_col]] == cl]
            if (length(pix) > 0) {
                all_intensities[[cl]] <- as.vector(itx_mat[, pix])
            }
        }
        df_all <- data.frame(
            intensity = unlist(all_intensities),
            cluster = rep(names(all_intensities), times = sapply(all_intensities, length)),
            stringsAsFactors = FALSE
        )
        df_all
    })
    
    violin_gg <- reactive({
        req(input$violin_mode)
        if (input$violin_mode == "single") {
            req(violin_df_single())
            df <- violin_df_single()
        } else {
            req(violin_df_all())
            df <- violin_df_all()
        }
        
        cluster_colors <- get_cluster_colors(df$cluster, pixel_clusters())
        
        if (input$violin_transpose) {
            p <- ggplot(df, aes(x = intensity, y = cluster, fill = cluster))
            x_label <- if(input$violin_transform == "log2") "Log2(intensity+1)" else "Raw intensity"
            y_label <- "Cluster"
            x_angle <- 0
            x_hjust <- 0.5
        } else {
            p <- ggplot(df, aes(x = cluster, y = intensity, fill = cluster))
            x_label <- "Cluster"
            y_label <- if(input$violin_transform == "log2") "Log2(intensity+1)" else "Raw intensity"
            x_angle <- input$violin_x_angle
            x_hjust <- 1
        }
        
        p <- p +
            scale_fill_manual(values = cluster_colors, drop = FALSE) +
            theme_bw() +
            labs(title = if(input$violin_mode == "single") input$violin_feature else "All features: global intensity distribution",
                 x = x_label, y = y_label) +
            theme(legend.position = "none",
                  axis.title = element_text(size = input$violin_axis_title_size),
                  axis.text.x = element_text(size = input$violin_axis_text_size, angle = x_angle, hjust = x_hjust),
                  axis.text.y = element_text(size = input$violin_axis_text_size),
                  plot.title = element_text(size = input$violin_title_size, hjust = 0.5))
        
        if (input$violin_plot_type == "violin") {
            p <- p + geom_violin(trim = TRUE)
        } else if (input$violin_plot_type == "boxplot") {
            p <- p + geom_boxplot(outlier.size = 0.3)
        } else if (input$violin_plot_type == "both") {
            p <- p + geom_violin(trim = TRUE, alpha = 0.6) + geom_boxplot(width = 0.1, outlier.size = 0.3)
        }
        
        p
    })
    
    output$violin_plot <- plotly::renderPlotly({
        req(violin_gg())
        plotly::ggplotly(violin_gg(), tooltip = c("x", "y", "fill")) %>%
            plotly::layout(dragmode = "lasso")
    })
    
    output$download_violin_png <- downloadHandler(
        filename = function() paste0("violin_plot_", Sys.Date(), ".png"),
        content = function(file) {
            ggsave(file, plot = violin_gg(), device = "png", width = 8, height = 6, dpi = 300)
        }
    )
    output$download_violin_pdf <- downloadHandler(
        filename = function() paste0("violin_plot_", Sys.Date(), ".pdf"),
        content = function(file) {
            ggsave(file, plot = violin_gg(), device = "pdf", width = 8, height = 6)
        }
    )
    
    observe({
        req(top_var())
        updateSelectizeInput(session, "violin_feature", choices = top_var(), server = TRUE)
    })
    
    observe({
        req(pixel_clusters())
        cluster_choices <- sort(unique(pixel_clusters()$hc_manual))
        updateSelectizeInput(session, "violin_clusters", choices = cluster_choices, server = TRUE)
    })
    
    # ----------------------------------------------------------------------------
    # PCA Plot (unchanged, uses helpers if needed)
    # ----------------------------------------------------------------------------
    pca_matrix <- eventReactive(input$run_pca, {
        withProgress(message = "PCA (fast irlba)", value = 0, {
            nbs_data <- neighbors()
            setProgress(0.05, detail = "Fetching data")
            if (!is.null(nbs_data) && !is.null(nbs_data[[2]])) {
                mat <- nbs_data[[2]]
                cluster_df <- nbs_data[[3]]
                cluster_df$cluster <- cluster_df$hc_orig
            } else {
                mat <- itx_res()
                cluster_df <- pixel_clusters()
                cluster_df$cluster <- cluster_df$hc_manual
            }
            if (!is.null(input$pca_clusters) && length(input$pca_clusters) > 0) {
                keep_pixels <- cluster_df$pixel_id[cluster_df$cluster %in% input$pca_clusters]
                mat <- mat[, intersect(keep_pixels, colnames(mat)), drop = FALSE]
                cluster_df <- cluster_df[cluster_df$pixel_id %in% colnames(mat), ]
            }
            if (ncol(mat) < 2) stop("Not enough pixels.")
            if (input$pca_transform == "log2") mat <- log2(mat + 1)
            if (inherits(mat, "DelayedArray")) mat <- as.matrix(mat)
            mat[!is.finite(mat)] <- 0
            
            setProgress(0.2, detail = "Removing constant features")
            sds <- matrixStats::rowSds(mat)
            keep <- which(sds > 0 & !is.na(sds))
            if (length(keep) < 2) stop("Not enough variable features")
            mat <- mat[keep, , drop = FALSE]
            
            setProgress(0.4, detail = "Selecting top variable features")
            if (nrow(mat) > input$pca_top_n) {
                sds <- matrixStats::rowSds(mat)
                top_idx <- order(sds, decreasing = TRUE)[1:input$pca_top_n]
                mat <- mat[top_idx, , drop = FALSE]
            }
            setProgress(0.6, detail = "Transposing & scaling")
            mat_t <- t(mat)
            mat_t[!is.finite(mat_t)] <- 0
            
            col_sds <- matrixStats::colSds(mat_t)
            const_cols <- which(col_sds == 0 | is.na(col_sds))
            if (length(const_cols) > 0) {
                mat_t <- mat_t[, -const_cols, drop = FALSE]
            }
            if (ncol(mat_t) < 2) stop("Too few variable pixels after removing constants.")
            
            setProgress(0.8, detail = "Running irlba")
            mat_scaled <- scale(mat_t, center = TRUE, scale = TRUE)
            mat_scaled[is.na(mat_scaled)] <- 0
            
            library(irlba)
            pca_res <- irlba(mat_scaled, nv = 10)
            scores <- pca_res$u %*% diag(pca_res$d)
            var_exp <- (pca_res$d^2) / sum(pca_res$d^2) * 100
            colnames(scores) <- paste0("PC", 1:ncol(scores))
            rownames(scores) <- rownames(mat_t)
            
            clusters <- cluster_df$cluster[match(rownames(mat_t), cluster_df$pixel_id)]
            setProgress(1, detail = "Done")
            list(pca = list(x = scores, sdev = pca_res$d), var_exp = var_exp, 
                 clusters = clusters, filter_clusters = input$pca_clusters)
        })
    })
    
    pca_gg <- reactive({
        req(pca_matrix())
        res <- pca_matrix()
        pc_x <- input$pca_dim_x
        pc_y <- input$pca_dim_y
        df_pca <- data.frame(PC1 = res$pca$x[, pc_x], PC2 = res$pca$x[, pc_y], Cluster = as.character(res$clusters))
        df_pca <- df_pca[!is.na(df_pca$Cluster), ]
        if (nrow(df_pca) == 0) return(NULL)
        
        # Use get_cluster_colors
        clusters_present <- unique(df_pca$Cluster)
        full_clusters <- pixel_clusters()$hc_manual
        full_colors <- pixel_clusters()$color
        color_map <- c()
        for (cl in clusters_present) {
            idx <- which(full_clusters == cl)
            if (length(idx) > 0) {
                tbl <- table(full_colors[idx])
                color_map[cl] <- names(tbl)[which.max(tbl)]
            }
        }
        missing <- setdiff(clusters_present, names(color_map))
        if (length(missing) > 0) {
            pal <- if (requireNamespace("khroma", quietly = TRUE)) khroma::color("smoothrainbow")(length(missing)) else rainbow(length(missing))
            names(pal) <- missing
            color_map <- c(color_map, pal)
        }
        
        filter_text <- if (!is.null(res$filter_clusters) && length(res$filter_clusters) > 0) {
            paste0(" (filtered to clusters: ", paste(res$filter_clusters, collapse=", "), ")")
        } else ""
        title_text <- paste0("Principal Component Analysis of Metabolite Intensities", filter_text)
        
        p <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Cluster)) +
            geom_point(size = input$pca_point_size, alpha = input$pca_point_alpha) +
            scale_color_manual(values = color_map) +
            theme_bw() +
            labs(x = paste0("PC", pc_x, " (", round(res$var_exp[pc_x], 1), "%)"),
                 y = paste0("PC", pc_y, " (", round(res$var_exp[pc_y], 1), "%)"),
                 title = title_text) +
            theme(axis.title = element_text(size = input$pca_axis_title_size),
                  axis.text = element_text(size = input$pca_axis_text_size),
                  plot.title = element_text(size = input$pca_title_size, hjust = 0.5),
                  legend.text = element_text(size = input$pca_legend_text_size),
                  legend.title = element_blank())
        if (input$pca_add_ellipse && length(unique(df_pca$Cluster)) > 1) {
            p <- p + stat_ellipse(level = 0.95, linetype = "dashed")
        }
        p
    })
    
    output$pca_plot <- plotly::renderPlotly({
        req(pca_gg())
        plotly::ggplotly(pca_gg(), tooltip = c("x", "y", "colour"))
    })
    
    output$pca_variance <- renderPrint({
        req(pca_matrix())
        var_exp <- pca_matrix()$var_exp
        cat("Variance explained by first 5 principal components:\n")
        for (i in 1:min(5, length(var_exp))) cat(sprintf("  %s: %.2f%%\n", paste0("PC", i), var_exp[i]))
    })
    
    output$download_pca_png <- downloadHandler(
        filename = function() paste0("pca_plot_", Sys.Date(), ".png"),
        content = function(file) ggsave(file, plot = pca_gg(), device = "png", width = 8, height = 6, dpi = 300)
    )
    output$download_pca_pdf <- downloadHandler(
        filename = function() paste0("pca_plot_", Sys.Date(), ".pdf"),
        content = function(file) ggsave(file, plot = pca_gg(), device = "pdf", width = 8, height = 6)
    )
    
    observe({
        req(pixel_clusters())
        cluster_choices <- sort(unique(pixel_clusters()$hc_manual))
        updateSelectizeInput(session, "pca_clusters", choices = cluster_choices, server = TRUE)
    })
    
    # ----------------------------------------------------------------------------
    # Scatter plot (uses get_cluster_colors helper)
    # ----------------------------------------------------------------------------
    scatter_data <- eventReactive(input$run_scatter, {
        req(input$scatter_x, input$scatter_y)
        
        withProgress(message = "Preparing scatter data", value = 0, {
            nbs_data <- neighbors()
            setProgress(0.1, detail = "Fetching data")
            
            if (!is.null(nbs_data) && !is.null(nbs_data[[2]])) {
                mat <- nbs_data[[2]]
                cluster_df <- nbs_data[[3]]
                cluster_df$cluster <- cluster_df$hc_orig
            } else {
                mat <- itx_res()
                cluster_df <- pixel_clusters()
                cluster_df$cluster <- cluster_df[["hc_manual"]]
            }
            
            if (!is.null(input$scatter_clusters) && length(input$scatter_clusters) > 0) {
                keep_pixels <- cluster_df$pixel_id[cluster_df$cluster %in% input$scatter_clusters]
                mat <- mat[, intersect(keep_pixels, colnames(mat)), drop = FALSE]
                cluster_df <- cluster_df[cluster_df$pixel_id %in% colnames(mat), ]
            }
            if (ncol(mat) < 2) stop("Not enough pixels after filtering.")
            
            if (input$scatter_transform == "log2") mat <- log2(mat + 1)
            if (inherits(mat, "DelayedArray")) mat <- as.matrix(mat)
            mat[!is.finite(mat)] <- 0
            
            x_vals <- mat[input$scatter_x, ]
            y_vals <- mat[input$scatter_y, ]
            valid <- !is.na(x_vals) & !is.na(y_vals)
            pixel_ids <- names(x_vals)[valid]
            x_vals <- x_vals[valid]
            y_vals <- y_vals[valid]
            
            clusters <- cluster_df$cluster[match(pixel_ids, cluster_df$pixel_id)]
            df <- data.frame(x = as.numeric(x_vals), y = as.numeric(y_vals), 
                             Cluster = as.character(clusters), stringsAsFactors = FALSE)
            df <- df[!is.na(df$Cluster) & df$Cluster != "", , drop = FALSE]
            
            setProgress(1, detail = "Done")
            list(df = df, x_name = input$scatter_x, y_name = input$scatter_y, filter_clusters = input$scatter_clusters)
        })
    })
    
    scatter_gg <- reactive({
        req(scatter_data())
        dat <- scatter_data()
        df <- dat$df
        if (nrow(df) == 0) {
            showNotification("No valid data points after filtering", type = "error")
            return(NULL)
        }
        
        cluster_colors <- get_cluster_colors(df$Cluster, pixel_clusters())
        filter_text <- if (!is.null(dat$filter_clusters) && length(dat$filter_clusters) > 0) {
            paste0(" (filtered to clusters: ", paste(dat$filter_clusters, collapse=", "), ")")
        } else ""
        title_text <- paste0("Scatter plot: ", dat$x_name, " vs ", dat$y_name, filter_text)
        
        p <- ggplot(df, aes(x = x, y = y, color = Cluster)) +
            geom_point(size = input$scatter_point_size, alpha = input$scatter_point_alpha) +
            scale_color_manual(values = cluster_colors, na.value = "grey50", drop = FALSE) +
            theme_bw() +
            labs(x = dat$x_name, y = dat$y_name, title = title_text) +
            theme(axis.title = element_text(size = input$scatter_axis_title_size),
                  axis.text = element_text(size = input$scatter_axis_text_size),
                  plot.title = element_text(size = input$scatter_title_size, hjust = 0.5),
                  legend.text = element_text(size = input$scatter_legend_text_size),
                  legend.title = element_blank())
        if (input$scatter_add_line) {
            p <- p + geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed")
        }
        p
    })
    
    output$scatter_plot <- plotly::renderPlotly({
        req(scatter_gg())
        plotly::ggplotly(scatter_gg(), tooltip = c("x", "y", "colour"))
    })
    
    output$download_scatter_png <- downloadHandler(
        filename = function() paste0("scatter_plot_", Sys.Date(), ".png"),
        content = function(file) ggsave(file, plot = scatter_gg(), device = "png", width = 8, height = 6, dpi = 300)
    )
    output$download_scatter_pdf <- downloadHandler(
        filename = function() paste0("scatter_plot_", Sys.Date(), ".pdf"),
        content = function(file) ggsave(file, plot = scatter_gg(), device = "pdf", width = 8, height = 6)
    )
    
    observe({
        req(top_var())
        updateSelectizeInput(session, "scatter_x", choices = top_var(), server = TRUE)
        updateSelectizeInput(session, "scatter_y", choices = top_var(), server = TRUE)
    })
    observe({
        req(pixel_clusters())
        cluster_choices <- sort(unique(pixel_clusters()$hc_manual))
        updateSelectizeInput(session, "scatter_clusters", choices = cluster_choices, server = TRUE)
    })
    
    # Populate cluster dropdown for Network tab (unchanged)
    observe({
        req(pixel_clusters())
        cluster_choices <- sort(unique(pixel_clusters()$hc_manual))
        updateSelectInput(session, "net_selected_domain",
                          choices = cluster_choices,
                          selected = if (length(cluster_choices) > 0) cluster_choices[1] else NULL)
    })
    
    # ----------------------------------------------------------------------------
    # UMAP Plot (unchanged)
    # ----------------------------------------------------------------------------
    umap_data <- eventReactive(input$run_umap, {
        withProgress(message = "UMAP (parallel via uwot)", value = 0, {
            nbs_data <- neighbors()
            setProgress(0.1, detail = "Preparing data")
            if (!is.null(nbs_data) && !is.null(nbs_data[[2]])) {
                mat <- nbs_data[[2]]
                cluster_df <- nbs_data[[3]]
                cluster_df$cluster <- cluster_df$hc_orig
            } else {
                mat <- itx_res()
                cluster_df <- pixel_clusters()
                cluster_df$cluster <- cluster_df$hc_manual
            }
            if (input$umap_transform == "log2") mat <- log2(mat + 1)
            if (inherits(mat, "DelayedArray")) mat <- as.matrix(mat)
            mat[!is.finite(mat)] <- 0
            setProgress(0.2, detail = "Removing constant features")
            sds <- matrixStats::rowSds(mat)
            keep <- which(sds > 0 & !is.na(sds))
            if (length(keep) < 2) stop("Not enough variable features")
            mat <- mat[keep, , drop = FALSE]
            setProgress(0.4, detail = "Selecting top features")
            if (nrow(mat) > input$umap_top_n) {
                sds <- matrixStats::rowSds(mat)
                top_idx <- order(sds, decreasing = TRUE)[1:input$umap_top_n]
                mat <- mat[top_idx, , drop = FALSE]
            }
            setProgress(0.6, detail = "Transposing")
            mat_t <- t(mat)
            mat_t[!is.finite(mat_t)] <- 0
            setProgress(0.7, detail = "Running UMAP (uwot)")
            set.seed(12345)
            umap_res <- uwot::umap(mat_t,
                                   n_neighbors = input$umap_n_neighbors,
                                   min_dist = input$umap_min_dist,
                                   n_components = max(input$umap_dim_x, input$umap_dim_y, 2),
                                   metric = "euclidean",
                                   verbose = FALSE,
                                   n_threads = availableCores() - 1)
            colnames(umap_res) <- paste0("UMAP", 1:ncol(umap_res))
            sample_names <- rownames(mat_t)
            clusters <- cluster_df$cluster[match(sample_names, cluster_df$pixel_id)]
            setProgress(1, detail = "Done")
            list(umap = umap_res, clusters = clusters)
        })
    })
    
    umap_gg <- reactive({
        req(umap_data())
        res <- umap_data()
        if (input$umap_swap_axes) {
            x_comp <- input$umap_dim_y
            y_comp <- input$umap_dim_x
        } else {
            x_comp <- input$umap_dim_x
            y_comp <- input$umap_dim_y
        }
        df_umap <- data.frame(
            UMAP1 = res$umap[, x_comp],
            UMAP2 = res$umap[, y_comp],
            Cluster = as.character(res$clusters)
        )
        df_umap <- df_umap[!is.na(df_umap$Cluster), ]
        if (nrow(df_umap) == 0) {
            showNotification("No valid clusters found for UMAP", type = "error")
            return(NULL)
        }
        clusters_present <- unique(df_umap$Cluster)
        full_clusters <- pixel_clusters()$hc_manual
        full_colors <- pixel_clusters()$color
        color_map <- c()
        for (cl in clusters_present) {
            idx <- which(full_clusters == cl)
            if (length(idx) > 0) {
                tbl <- table(full_colors[idx])
                color_map[cl] <- names(tbl)[which.max(tbl)]
            }
        }
        missing <- setdiff(clusters_present, names(color_map))
        if (length(missing) > 0) {
            pal <- if (requireNamespace("khroma", quietly = TRUE)) khroma::color("smoothrainbow")(length(missing)) else rainbow(length(missing))
            names(pal) <- missing
            color_map <- c(color_map, pal)
        }
        p <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
            geom_point(size = input$umap_point_size, alpha = input$umap_point_alpha) +
            scale_color_manual(values = color_map) +
            theme_bw() +
            labs(x = paste0("UMAP", x_comp), y = paste0("UMAP", y_comp),
                 title = "UMAP of Metabolite Intensities") +
            theme(axis.title = element_text(size = input$umap_axis_title_size),
                  axis.text = element_text(size = input$umap_axis_text_size),
                  plot.title = element_text(size = input$umap_title_size, hjust = 0.5),
                  legend.text = element_text(size = input$umap_legend_text_size),
                  legend.title = element_blank())
        if (input$umap_add_ellipse && length(unique(df_umap$Cluster)) > 1) {
            p <- p + stat_ellipse(level = 0.95, linetype = "dashed")
        }
        p
    })
    
    output$umap_plot <- plotly::renderPlotly({
        req(umap_gg())
        plotly::ggplotly(umap_gg(), tooltip = c("x", "y", "colour"))
    })
    
    output$umap_info <- renderPrint({
        req(umap_data())
        cat("UMAP parameters:\n")
        cat("  n_neighbors =", input$umap_n_neighbors, "\n")
        cat("  min_dist =", input$umap_min_dist, "\n")
        cat("  Input features =", ncol(umap_data()$umap), "samples\n")
    })
    
    output$download_umap_png <- downloadHandler(
        filename = function() paste0("umap_plot_", Sys.Date(), ".png"),
        content = function(file) ggsave(file, plot = umap_gg(), device = "png", width = 8, height = 6, dpi = 300)
    )
    output$download_umap_pdf <- downloadHandler(
        filename = function() paste0("umap_plot_", Sys.Date(), ".pdf"),
        content = function(file) ggsave(file, plot = umap_gg(), device = "pdf", width = 8, height = 6)
    )
    
    # ----------------------------------------------------------------------------
    # Spatial statistics (unchanged)
    # ----------------------------------------------------------------------------
    spatial_stats_res <- eventReactive(input$test_run_spatial_stats, {
        top_features <- as.integer(input$user_top_var)
        top_features <- top_var()[1:top_features]
        df_tmp <- data.frame(pixel_id = pixel_clusters()$pixel_id,
                             x_coord = pixel_clusters()$x_coord,
                             y_coord = pixel_clusters()$y_coord)
        itx_tmp <- itx_res()[top_features, df_tmp$pixel_id, drop = FALSE]
        
        withProgress(message = "Calculating spatial weights...", value = 0, {
            k <- if (nrow(df_tmp) >= 45000) 8 else 4
            spw <- calculate_spatial_weights(coords = df_tmp, k = k)
        })
        
        withProgress(message = "Calculating statistics in parallel (row-wise)", value = 0, {
            library(BiocParallel)
            BPPARAM <- MulticoreParam(workers = max(1, availableCores() - 1))
            
            if (inherits(itx_tmp, "DelayedArray")) {
                itx_tmp <- as.matrix(itx_tmp)
            }
            
            results <- bplapply(seq_len(nrow(itx_tmp)), function(i) {
                row_vals <- itx_tmp[i, , drop = FALSE]
                spatial_stats(x = row_vals, spw = spw)
            }, BPPARAM = BPPARAM)
            
            sp_res <- dplyr::bind_rows(results)
            sp_res
        })
    })
    
    output$sp_output <- DT::renderDT(
        DT::datatable(spatial_stats_res() %>%
                          dplyr::mutate(dplyr::across(c("moran_i", "geary_c"), round, digits = 5)) %>%
                          dplyr::rename(Molecule = feature, `Moran's I` = moran_i, `Geary's C` = geary_c),
                      options = list(scrollX = TRUE, scrollCollapse = TRUE, pageLength = 10))
    )
    
    # ----------------------------------------------------------------------------
    # Spatial gradients (unchanged)
    # ----------------------------------------------------------------------------
    spatial_gradients_res <- eventReactive(input$test_run_spatial_gradients, {
        top_features <- as.integer(input$user_top_var_grad)
        top_features <- top_var()[1:top_features]
        df_tmp <- data.frame(pixel_id = pixel_clusters()[['pixel_id']],
                             x_coord = pixel_clusters()[['x_coord']],
                             y_coord = pixel_clusters()[['y_coord']],
                             annots = pixel_clusters()[[input$annotation_test_gradients]])
        df_tmp[[4]] <- as.character(df_tmp[[4]])
        itx_tmp <- itx_res()[top_features, df_tmp[['pixel_id']]]
        min_nb <- 9
        
        withProgress(message = "Calculating spatial gradients...", value = 0.2, {
            spg_tests_df <- spatial_gradient(
                x = itx_tmp,
                sp_df = df_tmp[, 1:3],
                ref = df_tmp[['pixel_id']][df_tmp[[4]] == input$ref_group],
                distsumm = input$summ_type,
                nbs_ls = NULL,
                min_nb = min_nb,
                log_dist = TRUE
            )
            incProgress(0.8)
        })
        
        spg_tests_df
    })
    
    output$spgradient_output <- DT::renderDT({
        res <- spatial_gradients_res()
        if (is.null(res) || nrow(res) == 0 || !any(c("molecule", "lm_coef") %in% colnames(res))) {
            DT::datatable(
                data.frame(Message = "No spatial gradient results. Try different parameters."),
                options = list(scrollX = TRUE, dom = 't'),
                rownames = FALSE
            )
        } else {
            numeric_cols <- intersect(c("lm_coef", "lm_pval", "spearman_r", "spearman_r_pval", "spearman_r_pval_adj"), colnames(res))
            if (length(numeric_cols) > 0) {
                res <- res %>% dplyr::mutate(dplyr::across(dplyr::all_of(numeric_cols), round, digits = 5))
            }
            rename_map <- c(molecule = "Molecule", lm_coef = "LM Slope", lm_pval = "Slope p-value",
                            spearman_r = "Spearman R", spearman_r_pval = "Spearman p-value",
                            spearman_r_pval_adj = "FDR (Spearman)", pval_comment = "Comments")
            for (old in names(rename_map)) {
                if (old %in% colnames(res)) colnames(res)[colnames(res) == old] <- rename_map[old]
            }
            DT::datatable(res, options = list(scrollX = TRUE, scrollCollapse = TRUE, pageLength = 10))
        }
    })
    
    observe({
        req(input_validate(), pixel_clusters(), input$annotation_test_gradients)
        cat('EXTRACT LIST OF LABELS FOR USER SELECTION (GRADIENTS 2)\n')
        cluster_opts <- unique(pixel_clusters()[[input$annotation_test_gradients]])
        cluster_opts <- cluster_opts[!is.na(cluster_opts)]
        updateSelectizeInput(session, 'ref_group', 
                             choices = cluster_opts, 
                             server = TRUE, 
                             selected = if (length(cluster_opts) >= 3) cluster_opts[3] else cluster_opts[1])
    })
    
    # ========================= FPCA (uses helpers) =========================
    
    fpca_results <- eventReactive(input$run_fpca, {
        req(input_validate())
        
        nbs_data <- neighbors()
        use_collapsed <- input$fpca_use_collapsed && !is.null(nbs_data[[1]])
        if (use_collapsed) {
            cl_df <- nbs_data[[3]]
            coords <- cl_df[, c("x_coord", "y_coord")]
            cluster_assignments <- cl_df$hc_orig
            cat("Using collapsed neighbor data\n")
        } else {
            cl_df <- pixel_clusters()
            coords <- cl_df[, c("x_coord", "y_coord")]
            cluster_assignments <- cl_df$hc_manual
            cat("Using original pixel clusters\n")
        }
        
        subsample_rate <- input$fpca_subsample
        if (!is.null(subsample_rate) && subsample_rate > 1) {
            keep_idx <- seq(1, nrow(coords), by = subsample_rate)
            coords <- coords[keep_idx, ]
            cluster_assignments <- cluster_assignments[keep_idx]
            cat("Subsampled pixels: keep 1 of", subsample_rate, "→", nrow(coords), "pixels\n")
        }
        
        min_pixels <- 10
        valid_clusters <- names(which(table(cluster_assignments) >= min_pixels))
        if (length(valid_clusters) < 3) {
            showNotification("Need at least 3 clusters with ≥10 pixels each", type = "error")
            return(NULL)
        }
        
        if (!is.null(input$fpca_max_clusters) && input$fpca_max_clusters > 0) {
            cluster_sizes <- table(cluster_assignments)
            largest <- names(sort(cluster_sizes, decreasing = TRUE)[
                1:min(input$fpca_max_clusters, length(valid_clusters))
            ])
            valid_clusters <- intersect(valid_clusters, largest)
            cat("Limiting to", length(valid_clusters), "largest clusters\n")
        }
        
        r_seq <- seq(0, input$fpca_r_max, by = input$fpca_r_step)
        
        withProgress(message = "Computing spatial G‑functions", value = 0, {
            all_curves <- list()
            total <- length(valid_clusters)
            
            for (idx in seq_along(valid_clusters)) {
                cl <- valid_clusters[idx]
                incProgress(1/total, detail = paste("Cluster", cl))
                
                pix_idx <- which(cluster_assignments == cl)
                if (length(pix_idx) < min_pixels) next
                coords_cl <- coords[pix_idx, ]
                
                win <- spatstat.geom::owin(xrange = range(coords_cl$x_coord),
                                           yrange = range(coords_cl$y_coord))
                pp_cl <- spatstat.geom::ppp(x = coords_cl$x_coord, y = coords_cl$y_coord, window = win)
                
                obs <- tryCatch({
                    spatstat.explore::Gest(pp_cl, r = r_seq, rmax = input$fpca_r_max, correction = "rs")
                }, error = function(e) NULL)
                if (is.null(obs)) next
                obs_df <- data.frame(r = obs$r, observed = obs$rs)
                
                n_perms <- input$fpca_n_perms
                if (nrow(coords) > 50000) n_perms <- min(n_perms, 30)
                if (nrow(coords) > 100000) n_perms <- min(n_perms, 15)
                
                perm_curves <- lapply(seq_len(n_perms), function(i) {
                    set.seed(333 + i)
                    perm_pix <- sample(seq_len(nrow(coords)), size = length(pix_idx), replace = FALSE)
                    coords_perm <- coords[perm_pix, ]
                    win_perm <- spatstat.geom::owin(xrange = range(coords_perm$x_coord),
                                                    yrange = range(coords_perm$y_coord))
                    pp_perm <- spatstat.geom::ppp(x = coords_perm$x_coord, y = coords_perm$y_coord, window = win_perm)
                    g_perm <- spatstat.explore::Gest(pp_perm, r = r_seq, rmax = input$fpca_r_max, correction = "rs")
                    data.frame(r = g_perm$r, permed = g_perm$rs)
                }) %>% dplyr::bind_rows() %>%
                    dplyr::group_by(r) %>%
                    dplyr::summarise(permed = mean(permed, na.rm = TRUE), .groups = "drop")
                
                curve <- dplyr::left_join(obs_df, perm_curves, by = "r") %>%
                    dplyr::mutate(fundiff = observed - permed, cluster = cl) %>%
                    dplyr::filter(r > 0, is.finite(fundiff))
                all_curves[[cl]] <- curve
                gc()
            }
            
            if (length(all_curves) < 3) {
                showNotification("Not enough valid curves for FPCA", type = "error")
                return(NULL)
            }
            
            all_out_summ <- dplyr::bind_rows(all_curves)
            incProgress(0.8, detail = "Running FPCA")
            
            if ("hc_manual" %in% colnames(cl_df)) {
                metadata <- cl_df %>%
                    dplyr::select(cluster = hc_manual, color) %>%
                    dplyr::distinct() %>%
                    dplyr::mutate(cluster = as.character(cluster))
                spatial <- cl_df %>%
                    dplyr::select(cluster = hc_manual, X = x_coord, Y = y_coord) %>%
                    dplyr::mutate(cluster = as.character(cluster))
            } else {
                metadata <- cl_df %>%
                    dplyr::select(cluster = hc_orig, color) %>%
                    dplyr::distinct() %>%
                    dplyr::mutate(cluster = as.character(cluster))
                spatial <- cl_df %>%
                    dplyr::select(cluster = hc_orig, X = x_coord, Y = y_coord) %>%
                    dplyr::mutate(cluster = as.character(cluster))
            }
            metadata <- metadata %>% dplyr::filter(cluster %in% valid_clusters)
            spatial <- spatial %>% dplyr::filter(cluster %in% valid_clusters)
            
            obj <- make_mxfda(metadata = metadata,
                              spatial = spatial,
                              subject_key = "color",
                              sample_key = "cluster")
            
            obj <- add_summary_function(obj,
                                        summary_function_data = all_out_summ %>%
                                            dplyr::mutate(cluster = as.character(cluster)),
                                        metric = "uni g")
            
            obj <- run_fpca(obj,
                            metric = "uni g",
                            r = "r",
                            value = "fundiff",
                            pve = input$fpca_pve)
            
            # Use helper extract_fpca_safely
            fpca_list <- extract_fpca_safely(obj)
            if (is.null(fpca_list)) {
                showNotification("FPCA extraction failed – no scores found", type = "error")
                return(NULL)
            }
            
            if (is.null(rownames(fpca_list$scores)) && nrow(fpca_list$scores) == nrow(metadata)) {
                rownames(fpca_list$scores) <- metadata$cluster
            } else if (is.null(rownames(fpca_list$scores))) {
                rownames(fpca_list$scores) <- seq_len(nrow(fpca_list$scores))
            }
            
            list(
                scores = fpca_list$scores,
                varprop = fpca_list$varprop,
                metadata = metadata,
                fundiff_data = all_out_summ,
                clusters_used = valid_clusters
            )
        })
    })
    
    fpca_score_gg <- reactive({
        req(fpca_results())
        res <- fpca_results()
        if (is.null(res$scores) || ncol(res$scores) < 2) return(NULL)
        
        point_size <- ifnull(input$fpca_point_size, 3)
        label_size <- ifnull(input$fpca_label_size, 4)
        title_size <- ifnull(input$fpca_title_size, 16)
        axis_title_size <- ifnull(input$fpca_axis_title_size, 14)
        axis_text_size <- ifnull(input$fpca_axis_text_size, 12)
        legend_text_size <- ifnull(input$fpca_legend_text_size, 12)
        
        scores <- res$scores
        if (is.null(rownames(scores))) rownames(scores) <- seq_len(nrow(scores))
        colnames(scores) <- paste0("FPC", seq_len(ncol(scores)))
        varprop <- res$varprop
        
        scores_df <- data.frame(cluster = rownames(scores), scores, stringsAsFactors = FALSE)
        scores_df$cluster <- as.character(scores_df$cluster)
        res$metadata$cluster <- as.character(res$metadata$cluster)
        scores_df <- dplyr::left_join(scores_df, res$metadata, by = "cluster")
        
        if (!"color" %in% colnames(scores_df)) scores_df$color <- NA
        scores_df$color[is.na(scores_df$color)] <- "#CCCCCC"
        
        color_map <- tryCatch({
            get_cluster_colors(scores_df$cluster, pixel_clusters())
        }, error = function(e) {
            cl_uniq <- unique(scores_df$cluster)
            col_pal <- scales::hue_pal()(length(cl_uniq))
            names(col_pal) <- cl_uniq
            col_pal
        })
        
        p <- ggplot(scores_df, aes(x = FPC1, y = FPC2, color = cluster)) +
            geom_point(size = point_size, alpha = 0.8) +
            ggrepel::geom_text_repel(aes(label = cluster), show.legend = FALSE, size = label_size) +
            scale_color_manual(values = color_map) +
            theme_bw(base_size = axis_text_size) +
            labs(x = if (!is.null(varprop) && length(varprop) >= 1) paste0("FPC1 (", round(varprop[1]*100, 1), "%)") else "FPC1",
                 y = if (!is.null(varprop) && length(varprop) >= 2) paste0("FPC2 (", round(varprop[2]*100, 1), "%)") else "FPC2",
                 title = "FPCA of spatial G‑function differences",
                 color = "Cluster") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
                  axis.title = element_text(size = axis_title_size),
                  axis.text = element_text(size = axis_text_size),
                  legend.title = element_blank(),
                  legend.text = element_text(size = legend_text_size),
                  legend.position = "right")
        p
    })
    
    fpca_var_gg <- reactive({
        req(fpca_results())
        varprop <- fpca_results()$varprop
        if (is.null(varprop) || length(varprop) == 0) return(NULL)
        
        title_size <- ifnull(input$fpca_var_title_size, 14)
        axis_title_size <- ifnull(input$fpca_var_axis_title_size, 12)
        axis_text_size <- ifnull(input$fpca_var_axis_text_size, 10)
        bar_width <- ifnull(input$fpca_var_bar_width, 0.7)
        line_type <- ifnull(input$fpca_line_type, "individual")
        secondary <- ifnull(input$fpca_var_secondary, FALSE)
        
        n_pcs <- length(varprop)
        var_df <- data.frame(
            PC = factor(1:n_pcs, levels = 1:n_pcs),
            Variance = varprop * 100,
            Cumulative = cumsum(varprop) * 100
        )
        
        p <- ggplot(var_df, aes(x = PC)) +
            geom_col(aes(y = Variance), fill = "steelblue", alpha = 0.7, width = bar_width) +
            labs(x = "Principal Component", y = "Variance explained (%)",
                 title = "Variance explained by FPCs") +
            theme_bw(base_size = axis_text_size) +
            theme(plot.title = element_text(hjust = 0.5, size = title_size),
                  axis.title = element_text(size = axis_title_size),
                  axis.text = element_text(size = axis_text_size),
                  panel.grid.minor = element_blank())
        
        if (line_type == "individual") {
            p <- p + geom_line(aes(y = Variance, group = 1), color = "red", size = 1, linetype = "solid") +
                geom_point(aes(y = Variance), color = "red", size = 2) +
                labs(caption = "Red line: variance explained by each PC")
        } else {
            if (secondary) {
                p <- p + geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
                    geom_point(aes(y = Cumulative), color = "red", size = 2) +
                    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Cumulative variance (%)"))
            } else {
                p <- p + geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1, linetype = "dashed") +
                    geom_point(aes(y = Cumulative), color = "red", size = 2) +
                    labs(caption = "Red dashed line: cumulative variance")
            }
        }
        p
    })
    
    output$fpca_score_plot <- plotly::renderPlotly({
        req(fpca_score_gg())
        plotly::ggplotly(fpca_score_gg(), tooltip = c("x", "y", "colour", "label")) %>%
            plotly::layout(dragmode = "lasso")
    })
    
    output$fpca_var_plot <- plotly::renderPlotly({
        req(fpca_var_gg())
        plotly::ggplotly(fpca_var_gg(), tooltip = c("x", "y")) %>%
            plotly::layout(dragmode = "lasso")
    })
    
    output$download_fpca_score_png <- downloadHandler(
        filename = function() paste0("fpca_score_plot_", Sys.Date(), ".png"),
        content = function(file) ggsave(file, plot = fpca_score_gg(), device = "png", width = 8, height = 6, dpi = 300)
    )
    output$download_fpca_score_pdf <- downloadHandler(
        filename = function() paste0("fpca_score_plot_", Sys.Date(), ".pdf"),
        content = function(file) ggsave(file, plot = fpca_score_gg(), device = "pdf", width = 8, height = 6)
    )
    
    output$download_fpca_var_png <- downloadHandler(
        filename = function() paste0("fpca_variance_plot_", Sys.Date(), ".png"),
        content = function(file) ggsave(file, plot = fpca_var_gg(), device = "png", width = 8, height = 6, dpi = 300)
    )
    output$download_fpca_var_pdf <- downloadHandler(
        filename = function() paste0("fpca_variance_plot_", Sys.Date(), ".pdf"),
        content = function(file) ggsave(file, plot = fpca_var_gg(), device = "pdf", width = 8, height = 6)
    )
    
    output$fpca_summary <- renderPrint({
        req(fpca_results())
        res <- fpca_results()
        if (is.null(res$scores)) {
            cat("FPCA did not produce valid results.\n")
            return()
        }
        cat("Number of FPCs retained:", ncol(res$scores), "\n")
        cat("Total variance explained:", round(sum(res$varprop) * 100, 1), "%\n")
        cat("Number of clusters used:", length(res$clusters_used), "\n")
    })
    
    output$download_fpca_scores <- downloadHandler(
        filename = function() paste0("fpca_scores_", Sys.Date(), ".csv"),
        content = function(file) {
            req(fpca_results())
            res <- fpca_results()
            if (is.null(res$scores)) return()
            scores <- res$scores
            if (is.null(rownames(scores))) rownames(scores) <- seq_len(nrow(scores))
            colnames(scores) <- paste0("FPC", seq_len(ncol(scores)))
            out <- data.frame(cluster = rownames(scores), scores, stringsAsFactors = FALSE)
            out <- dplyr::left_join(out, res$metadata, by = "cluster")
            write.csv(out, file, row.names = FALSE)
        }
    )
    
    # ----------------------------------------------------------------------------
    # METABOLITE NETWORK (uses plotly_network helper)
    # ----------------------------------------------------------------------------
    current_layout <- reactiveVal(NULL)
    
    network_data <- eventReactive(input$run_network, {
        req(itx_res(), pixel_clusters(), top_var())
        
        withProgress(message = "Building correlation network...", value = 0, {
            
            if (input$net_filter_mode == "global") {
                top_n <- input$net_top_n
                features <- top_var()[1:min(top_n, length(top_var()))]
            } else {
                dom <- input$net_selected_domain
                req(dom, pixel_clusters(), itx_res())
                cl_df <- pixel_clusters()
                dom_pixels <- cl_df$pixel_id[cl_df$hc_manual == dom]
                if (length(dom_pixels) < 3) {
                    showNotification("Selected domain has too few pixels (<3)", type = "warning")
                    return(NULL)
                }
                mat_dom <- itx_res()[, dom_pixels, drop = FALSE]
                if (inherits(mat_dom, "DelayedArray")) mat_dom <- as.matrix(mat_dom)
                row_var <- apply(mat_dom, 1, var, na.rm = TRUE)
                row_var <- row_var[is.finite(row_var) & row_var > 0]
                top_features <- names(sort(row_var, decreasing = TRUE))
                top_n <- input$net_top_per_domain
                features <- top_features[1:min(top_n, length(top_features))]
            }
            
            req(length(features) >= 3)
            incProgress(0.1, detail = paste("Selected", length(features), "features"))
            
            mat <- itx_res()[features, , drop = FALSE]
            if (inherits(mat, "DelayedArray")) mat <- as.matrix(mat)
            mat[!is.finite(mat)] <- 0
            
            if (ncol(mat) > 5000) {
                set.seed(123)
                samp <- sample(ncol(mat), min(5000, ncol(mat)))
                mat <- mat[, samp]
                incProgress(0.2, detail = "Subsampled to 5000 pixels")
            }
            
            incProgress(0.3, detail = "Computing correlation matrix...")
            cor_mat <- WGCNA::cor(t(mat), use = "pairwise.complete.obs", method = input$net_corr_method)
            WGCNA::enableWGCNAThreads(nThreads = availableCores() - 1)
            cor_mat[is.na(cor_mat)] <- 0
            
            incProgress(0.6, detail = "Building graph and filtering...")
            library(igraph)
            g <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", 
                                             weighted = TRUE, diag = FALSE)
            threshold <- input$net_corr_threshold
            g <- delete_edges(g, which(abs(E(g)$weight) < threshold))
            g <- delete_vertices(g, which(degree(g) < input$net_min_degree))
            
            if (vcount(g) < 3) {
                showNotification("Not enough nodes after filtering. Lower threshold or include more features.", 
                                 type = "warning", duration = 10)
                return(NULL)
            }
            
            incProgress(0.8, detail = "Preparing layout...")
            layout_choice <- switch(input$net_layout,
                                    fr = layout_with_fr(g),
                                    kk = layout_with_kk(g),
                                    circle = layout_in_circle(g))
            coords <- as.data.frame(layout_choice)
            colnames(coords) <- c("x", "y")
            
            nodes <- data.frame(
                id = V(g)$name,
                label = V(g)$name,
                degree = degree(g),
                x = coords$x,
                y = coords$y,
                stringsAsFactors = FALSE
            )
            
            edges <- data.frame(
                from = get.edgelist(g)[,1],
                to = get.edgelist(g)[,2],
                weight = E(g)$weight,
                stringsAsFactors = FALSE
            )
            
            incProgress(1, detail = "Done")
            
            current_layout(NULL)
            showNotification("New network built. Use 'Save Current Layout' after dragging nodes.", 
                             type = "message")
            
            list(nodes = nodes, edges = edges, graph = g)
        })
    })
    
    observeEvent(input$save_layout, {
        req(network_data())
        req(input$net_plot_type == "vis")
        visNetworkProxy("network_vis") %>%
            visGetPositions(input = "network_vis_positions")
        showNotification("Fetching current layout...", type = "message", duration = 2)
    })
    
    observeEvent(input$network_vis_positions, {
        pos <- input$network_vis_positions
        if (!is.null(pos) && length(pos) > 0) {
            coords <- data.frame(
                id = names(pos),
                x = sapply(pos, function(p) p$x),
                y = sapply(pos, function(p) p$y),
                stringsAsFactors = FALSE
            )
            dat <- network_data()
            if (!is.null(dat)) {
                nodes <- dat$nodes
                nodes <- merge(nodes, coords, by = "id", all.x = TRUE, suffixes = c("", ".new"))
                nodes$x <- ifelse(is.na(nodes$x.new), nodes$x, nodes$x.new)
                nodes$y <- ifelse(is.na(nodes$y.new), nodes$y, nodes$y.new)
                nodes <- nodes[, c("id", "x", "y")]
                current_layout(nodes)
                showNotification("Layout saved! Downloads will now use this arrangement.", 
                                 type = "message")
            }
        } else {
            showNotification("No positions received. Please drag nodes first, then save.", 
                             type = "warning")
        }
    })
    
    output$network_plot <- renderUI({
        req(network_data())
        if (input$net_plot_type == "vis") {
            visNetwork::visNetworkOutput("network_vis", height = "750px")
        } else {
            plotly::plotlyOutput("network_plotly", height = "750px")
        }
    })
    
    output$network_vis <- renderVisNetwork({
        req(network_data())
        dat <- network_data()
        if (is.null(dat)) {
            return(visNetwork(
                data.frame(id = 1, label = "No network generated"),
                data.frame(from = numeric(0), to = numeric(0))
            ))
        }
        
        nodes <- data.frame(
            id = dat$nodes$id,
            label = dat$nodes$label,
            size = dat$nodes$degree * 2 + 5,
            color = input$node_color,
            stringsAsFactors = FALSE
        )
        
        edges <- data.frame(
            from = dat$edges$from,
            to = dat$edges$to,
            weight = dat$edges$weight,
            color = ifelse(dat$edges$weight > 0, input$edge_color_pos, input$edge_color_neg),
            width = abs(dat$edges$weight) * 3,
            stringsAsFactors = FALSE
        )
        if (input$net_show_edge_labels) {
            edges$label <- as.character(round(edges$weight, 3))
        } else {
            edges$label <- NA_character_
        }
        
        layout_map <- switch(input$net_layout,
                             fr = "layout_with_fr",
                             kk = "layout_with_kk",
                             circle = "layout_in_circle")
        
        visNetwork(nodes, edges) %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
            visPhysics(enabled = input$net_interactive) %>%
            visIgraphLayout(layout = layout_map) %>%
            visEdges(arrows = "",
                     font = list(size = input$net_edge_label_size, color = "black"),
                     label = edges$label) %>%
            visNodes(font = list(size = input$net_node_label_size)) %>%
            visLayout(randomSeed = 123)
    })
    
    output$network_plotly <- renderPlotly({
        req(network_data())
        dat <- network_data()
        if (is.null(dat)) {
            return(plotly::plot_ly() %>%
                       layout(annotations = list(text = "No network generated")))
        }
        
        # Use helper plotly_network
        plotly_network(
            dat,
            show_labels = input$net_show_edge_labels,
            label_size = input$net_edge_label_size,
            node_label_size = input$net_node_label_size,
            node_color = input$node_color,
            edge_pos_color = input$edge_color_pos,
            edge_neg_color = input$edge_color_neg
        )
    })
    
    output$network_summary <- renderPrint({
        req(network_data())
        dat <- network_data()
        if (is.null(dat)) {
            cat("No network generated.\nTry adjusting thresholds or selecting more features.")
            return()
        }
        g <- dat$graph
        cat("Network summary:\n")
        cat("  Nodes:", vcount(g), "\n")
        cat("  Edges:", ecount(g), "\n")
        cat("  Average degree:", round(mean(degree(g)), 2), "\n")
        cat("  Density:", round(edge_density(g), 4), "\n")
        cat("  Positive edges:", sum(E(g)$weight > 0), "\n")
        cat("  Negative edges:", sum(E(g)$weight < 0), "\n")
    })
    
    observeEvent(input$screenshot_btn, {
        plot_id <- if (input$net_plot_type == "vis") "network_vis" else "network_plotly"
        js_code <- sprintf("
    html2canvas(document.getElementById('%s'), {
      scale: 2,
      useCORS: true,
      allowTaint: true,
      backgroundColor: '#ffffff'
    }).then(function(canvas) {
      var link = document.createElement('a');
      link.download = 'network_screenshot.png';
      link.href = canvas.toDataURL('image/png');
      link.click();
    }).catch(function(error) {
      Shiny.setInputValue('screenshot_error', error.message);
    });
  ", plot_id)
        shinyjs::runjs(js_code)
    })
    
    observeEvent(input$screenshot_error, {
        showNotification(paste("Screenshot failed:", input$screenshot_error), type = "error")
    })
    
    output$download_network_png <- downloadHandler(
        filename = function() paste0("metabolite_network_", Sys.Date(), ".png"),
        content = function(file) {
            req(network_data())
            dat <- network_data()
            if (is.null(dat)) return()
            library(igraph)
            g <- dat$graph
            V(g)$color <- input$node_color
            E(g)$color <- ifelse(E(g)$weight > 0, input$edge_color_pos, input$edge_color_neg)
            E(g)$width <- abs(E(g)$weight) * 3
            
            if (!is.null(current_layout())) {
                layout_mat <- as.matrix(current_layout()[, c("x", "y")])
            } else {
                layout_mat <- as.matrix(dat$nodes[, c("x", "y")])
            }
            
            png(file, width = 1200, height = 1000, res = 120)
            plot(g,
                 vertex.color = V(g)$color,
                 vertex.size = degree(g) * 2 + 5,
                 edge.color = E(g)$color,
                 edge.width = E(g)$width,
                 edge.label = if (input$net_show_edge_labels) round(E(g)$weight, 3) else NA,
                 edge.label.cex = 0.6,
                 layout = layout_mat,
                 main = "Metabolite Correlation Network")
            dev.off()
        }
    )
    
    output$download_network_pdf <- downloadHandler(
        filename = function() paste0("metabolite_network_", Sys.Date(), ".pdf"),
        content = function(file) {
            req(network_data())
            dat <- network_data()
            if (is.null(dat)) return()
            library(igraph)
            g <- dat$graph
            V(g)$color <- input$node_color
            E(g)$color <- ifelse(E(g)$weight > 0, input$edge_color_pos, input$edge_color_neg)
            E(g)$width <- abs(E(g)$weight) * 3
            
            if (!is.null(current_layout())) {
                layout_mat <- as.matrix(current_layout()[, c("x", "y")])
            } else {
                layout_mat <- as.matrix(dat$nodes[, c("x", "y")])
            }
            
            pdf(file, width = 12, height = 10)
            plot(g,
                 vertex.color = V(g)$color,
                 vertex.size = degree(g) * 2 + 5,
                 edge.color = E(g)$color,
                 edge.width = E(g)$width,
                 edge.label = if (input$net_show_edge_labels) round(E(g)$weight, 3) else NA,
                 edge.label.cex = 0.6,
                 layout = layout_mat,
                 main = "Metabolite Correlation Network")
            dev.off()
        }
    )
    
    # ----------------------------------------------------------------------------
    # Domain Preferential Metabolites (uses helpers)
    # ----------------------------------------------------------------------------
    
    test_results <- eventReactive(input$run_occurrence, {
        req(input_validate(), pixel_clusters(), itx_res(), top_var())
        
        withProgress(message = "Computing tests...", value = 0, {
            n_features        <- min(as.integer(input$occ_top_n), length(top_var()))
            selected_features <- top_var()[1:n_features]
            
            cl_df   <- pixel_clusters()
            itx_mat <- itx_res()
            
            common_pixels <- intersect(cl_df$pixel_id, colnames(itx_mat))
            cl_df   <- cl_df[cl_df$pixel_id %in% common_pixels, ]
            itx_sub <- itx_mat[selected_features, cl_df$pixel_id, drop = FALSE]
            if (inherits(itx_sub, "DelayedArray")) itx_sub <- as.matrix(itx_sub)
            
            domains       <- sort(unique(cl_df$hc_manual))
            domain_labels <- cl_df$hc_manual
            test_method   <- input$occ_test_method
            incProgress(0.3, detail = paste("Testing method:", test_method))
            
            occ_list <- list()
            
            if (test_method == "wilcox") {
                for (dom in domains) {
                    in_dom <- domain_labels == dom
                    if (sum(in_dom) < 3 || sum(!in_dom) < 3) next
                    
                    p_vals <- apply(itx_sub, 1, function(x) {
                        wilcox.test(x[in_dom], x[!in_dom], exact = FALSE)$p.value
                    })
                    p_adj <- p.adjust(p_vals, method = "fdr")
                    log2fc <- compute_log2fc(itx_sub, domain_labels, dom)   # helper
                    
                    occ_list[[dom]] <- data.frame(
                        domain = dom,
                        metabolite = rownames(itx_sub),
                        log2fc = log2fc,
                        raw_p = p_vals,
                        adj_p = p_adj,
                        stringsAsFactors = FALSE
                    )
                }
            } else {  # permutation
                n_perm <- as.integer(input$occ_n_perm)
                incProgress(0.5, detail = paste("Permutations:", n_perm))
                
                for (dom in domains) {
                    in_dom <- domain_labels == dom
                    if (sum(in_dom) < 2 || sum(!in_dom) < 2) next
                    
                    obs_fc <- compute_log2fc(itx_sub, domain_labels, dom)
                    
                    perm_fc_mat <- vapply(seq_len(n_perm), function(i) {
                        compute_log2fc(itx_sub, sample(domain_labels), dom)
                    }, numeric(length(selected_features)))
                    
                    p_vals <- rowMeans(perm_fc_mat >= obs_fc, na.rm = TRUE)
                    p_vals <- (p_vals * n_perm + 1) / (n_perm + 1)
                    p_adj <- p.adjust(p_vals, method = "fdr")
                    
                    occ_list[[dom]] <- data.frame(
                        domain = dom,
                        metabolite = rownames(itx_sub),
                        log2fc = obs_fc,
                        raw_p = p_vals,
                        adj_p = p_adj,
                        stringsAsFactors = FALSE
                    )
                }
            }
            
            incProgress(0.9, detail = "Organizing results")
            results <- dplyr::bind_rows(occ_list)
            
            color_map <- get_cluster_colors(unique(results$domain), pixel_clusters())
            results$color <- color_map[as.character(results$domain)]
            results$color[is.na(results$color)] <- "#CCCCCC"
            
            incProgress(1, detail = "Done")
            results
        })
    })
    
    preferential_data <- reactive({
        req(test_results())
        
        res <- test_results()
        log2fc_thresh <- input$occ_log2fc_thresh
        fdr_thresh    <- if (input$occ_test_method == "wilcox") input$occ_fdr_thresh else input$occ_pval_thresh
        
        enriched <- !is.na(res$adj_p) & !is.na(res$log2fc) &
            abs(res$log2fc) >= log2fc_thresh & res$adj_p <= fdr_thresh
        
        res$enriched <- enriched
        
        occ_df <- res %>%
            dplyr::group_by(domain) %>%
            dplyr::summarise(
                n_metabolites = sum(enriched, na.rm = TRUE),
                n_pixels = NA_integer_,
                median_log2fc = median(log2fc[is.finite(log2fc)], na.rm = TRUE),
                min_raw_p = min(raw_p, na.rm = TRUE),
                min_adj_p = min(adj_p, na.rm = TRUE),
                preferential_metabolites = paste(metabolite[enriched], collapse = "; "),
                color = dplyr::first(color),
                .groups = "drop"
            ) %>%
            dplyr::mutate(
                preferential_metabolites = ifelse(n_metabolites == 0, "", preferential_metabolites)
            )
        
        cl_df <- pixel_clusters()
        pixel_counts <- table(cl_df$hc_manual)
        occ_df$n_pixels <- as.integer(pixel_counts[as.character(occ_df$domain)])
        occ_df$n_pixels[is.na(occ_df$n_pixels)] <- 0
        
        occ_df$sig_label <- dplyr::case_when(
            occ_df$n_metabolites == 0                         ~ NA_character_,
            is.na(occ_df$min_adj_p)                           ~ NA_character_,
            occ_df$min_adj_p < 0.001                          ~ "***",
            occ_df$min_adj_p < 0.01                           ~ "**",
            occ_df$min_adj_p < fdr_thresh                     ~ "*",
            TRUE                                              ~ NA_character_
        )
        
        occ_df
    })
    
    occ_ggplot <- reactive({
        show_stars <- input$occ_show_signif
        req(preferential_data())
        
        df <- preferential_data()
        if (nrow(df) == 0) {
            return(ggplot() + 
                       annotate("text", x=0.5, y=0.5, label="No metabolites meet criteria") + 
                       theme_void())
        }
        
        if (isTRUE(input$occ_sort_clusters)) {
            df$domain <- factor(df$domain, levels = df$domain[order(df$n_metabolites, decreasing = TRUE)])
        } else {
            df$domain <- factor(df$domain, levels = sort(unique(df$domain)))
        }
        
        color_vals <- setNames(df$color, as.character(df$domain))
        y_max   <- max(df$n_metabolites, na.rm = TRUE)
        y_nudge <- max(1, y_max * 0.06)
        
        fdr_thresh <- if (input$occ_test_method == "wilcox") input$occ_fdr_thresh else input$occ_pval_thresh
        test_label <- paste0("FDR ≤ ", fdr_thresh)
        
        p <- ggplot(df, aes(x = domain, y = n_metabolites, fill = domain)) +
            geom_col(width = 0.7, alpha = input$occ_bar_alpha) +
            scale_fill_manual(values = color_vals) +
            labs(x = "Domain / Cluster",
                 y = paste0("Log2FC > ", input$occ_log2fc_thresh, ", ", test_label),
                 title = "Domain preferential metabolites") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.title = element_text(size = input$occ_axis_title_size),
                  axis.text.x = element_text(size = input$occ_axis_text_size,
                                             angle = input$occ_x_angle,
                                             hjust = ifelse(input$occ_x_angle > 0, 1, 0.5)),
                  axis.text.y = element_text(size = input$occ_axis_text_size),
                  plot.title = element_text(size = input$occ_title_size, hjust = 0.5))
        
        if (isTRUE(input$occ_show_labels)) {
            p <- p + geom_text(aes(label = n_metabolites), vjust = -0.4, size = input$occ_axis_text_size / 3)
        }
        
        if (show_stars) {
            df_stars <- df[!is.na(df$sig_label) & df$n_metabolites > 0, ]
            if (nrow(df_stars) > 0) {
                p <- p + geom_text(
                    data = df_stars,
                    aes(y = n_metabolites + y_nudge, label = sig_label, color = sig_label == "ns"),
                    size = input$occ_axis_text_size / 2.5,
                    fontface = "bold",
                    show.legend = FALSE
                ) +
                    scale_color_manual(values = c("TRUE" = "grey60", "FALSE" = "black"))
            }
        }
        
        p
    })
    
    output$occ_plot <- plotly::renderPlotly({
        req(occ_ggplot())
        ggplotly(occ_ggplot(), tooltip = "text") %>%
            plotly::layout(hoverlabel = list(bgcolor = "white", font = list(size = 12)),
                           margin = list(b = 80))
    })
    
    output$occ_table <- DT::renderDT({
        req(preferential_data())
        df <- preferential_data()
        df <- df[, c("domain", "n_pixels", "n_metabolites", "median_log2fc", 
                     "min_raw_p", "min_adj_p", "preferential_metabolites")]
        DT::datatable(df, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    })
    
    output$download_occ_csv <- downloadHandler(
        filename = function() paste0("domain_preferential_metabolites_", Sys.Date(), ".csv"),
        content = function(file) {
            write.csv(preferential_data(), file, row.names = FALSE)
        }
    )
    
    output$download_occ_pdf <- downloadHandler(
        filename = function() paste0("domain_preferential_metabolites_", Sys.Date(), ".pdf"),
        content = function(file) {
            ggsave(file, plot = occ_ggplot(), width = 10, height = 7, dpi = 300)
        }
    )
    
    output$download_occ_png <- downloadHandler(
        filename = function() paste0("domain_preferential_metabolites_", Sys.Date(), ".png"),
        content = function(file) {
            ggsave(file, plot = occ_ggplot(), width = 10, height = 7, dpi = 300)
        }
    )
    
}) # END SERVER