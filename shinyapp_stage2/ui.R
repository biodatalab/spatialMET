library(shiny)
library(shinydashboard)
library(plotly)
library(colourpicker)
library(shinyjs) 

options(shiny.maxRequestSize = 10000 * 1024^2)

dashboardPage(
    dashboardHeader(title = "spatialMET"),
    dashboardSidebar(
        sidebarMenu(
            # ---- Data Import (single menu with two subitems) ----
            menuItem("Import data", 
                     icon = icon("upload"),
                     startExpanded = TRUE,
                     menuSubItem("Upload preprocessed data", tabName = "import", icon = icon("file-upload")),
                     menuSubItem("Preprocess raw data", tabName = "preprocess", icon = icon("cogs"))
            ),
            
            # ---- Analysis tools (shown only when data is loaded) ----
            conditionalPanel(
                condition = "output.fileUploaded || output.dataLoaded",
                sidebarMenu(
                    menuItem("Spatial Visualization",
                             icon = icon("map-marker"),
                             startExpanded = TRUE,
                             menuSubItem("Spatial domains", tabName = "spatialdomains", icon = icon("braille")),
                             menuSubItem("Intensities", tabName = "intensities", icon = icon("bolt"))
                    ),
                    hr(),
                    menuItem("Univariate Analysis",
                             icon = icon("chart-line"),
                             startExpanded = TRUE,
                             menuSubItem("Violin plot", tabName = "violin", icon = icon("chart-line")),
                             menuSubItem("PCA", tabName = "pca", icon = icon("chart-line")),
                             menuSubItem("Scatter plot", tabName = "scatter", icon = icon("chart-line")),
                             menuSubItem("UMAP", tabName = "umap", icon = icon("chart-line"))
                    ),
                    hr(),
                    menuItem("Spatial Statistics & Gradients",
                             icon = icon("map"),
                             startExpanded = TRUE,
                             menuSubItem("Spatial statistics", tabName = "spatial_stats", icon = icon("map")),
                             menuSubItem("Spatial gradient tests", tabName = "spatial_gradients", icon = icon("stairs")),
                             menuSubItem("FPCA (Spatial G‑function)", tabName = "fpca", icon = icon("chart-line"))
                    ),
                    hr(),
                    menuItem("Differential Abundance",
                             icon = icon("chart-simple"),
                             tabName = "diff_abund"),
                    hr(),
                    menuItem("Network & Enrichment",
                             icon = icon("project-diagram"),
                             startExpanded = TRUE,
                             menuSubItem("Metabolite Network", tabName = "met_network", icon = icon("project-diagram")),
                             menuSubItem("Metabolite Occurrence", tabName = "met_occurrence", icon = icon("chart-bar"))
                    )
                )
            )
        )
    ),
    dashboardBody(
        useShinyjs(),
        tags$head(
            tags$script(src = "https://cdn.jsdelivr.net/npm/html2canvas@1.4.1/dist/html2canvas.min.js")
        ),
        tabItems(
            # ========================= IMPORT DATA (upload preprocessed) =========================
            tabItem(
                tabName = "import",
                h1("spatialMET: Visualization and analysis of spatial metabolomics experiments", align = "center"),
                box(
                    width = 9,
                    style = "height: calc(100vh - 180px); overflow-y: auto;",
                    uiOutput("getting_started")
                ),
                box(
                    width = 3,
                    fileInput("hcdist_user_file", "Upload preprocessed data",
                              multiple = FALSE, accept = ".txt"),
                    fileInput("itx_user_file", "Upload a intensity table file",
                              multiple = FALSE, accept = ".tsv"),
                    textInput(inputId = "collapse_par", label = "", placeholder = NULL, value = "", width = "0px"),
                    fileInput("img_user_file", "Upload a tissue image (optional)",
                              multiple = FALSE, accept = "image/*")
                )
            ),
            
            # ========================= PREPROCESS RAW DATA =========================
            tabItem(
                tabName = "preprocess",
                h1("Preprocess raw imzML data", align = "center"),
                fluidRow(
                    box(
                        width = 4,
                        fileInput("raw_imzml_files", 
                                  "Upload .imzML and .ibd files (or a ZIP archive)",
                                  accept = c(".imzML", ".imzml", ".ibd", ".zip"),
                                  multiple = TRUE),
                        helpText("If uploading separately, select both .imzML and .ibd files together. 
            Alternatively, upload a ZIP archive containing both."),
                        hr(),
                        h4("Basic parameters"),
                        numericInput("prep_cores", "Number of cores", 
                                     value = max(1, parallel::detectCores() - 1), min = 1, max = 64, step = 1),
                        # REMOVED prep_n_clusters
                        numericInput("prep_n_mz_clusters", "Number of metabolite clusters",
                                     value = 50, min = 10, max = 200, step = 5),
                        hr(),
                        h4("Advanced options (optional)"),
                        tags$details(
                            tags$summary("Show filtering and clustering parameters"),
                            fluidRow(
                                column(6,
                                       # ------------------ New filtering inputs ------------------
                                       numericInput("prep_filter_present", "Filter: present fraction", 
                                                    value = 0.01, min = 0, max = 1, step = 0.01),
                                       numericInput("prep_filter_unlog_sd", "Filter: raw SD threshold", 
                                                    value = 250, min = 0, step = 10),
                                       numericInput("prep_filter_log2_sd", "Filter: log2 SD threshold", 
                                                    value = 0.5, min = 0, step = 0.1),
                                       # ----------------------------------------------------------
                                       numericInput("prep_minkowski_p", "Minkowski distance p", 
                                                    value = 1.5, min = 1, max = 5, step = 0.1),
                                       numericInput("prep_distpow", "Distance power", 
                                                    value = 1.5, min = 1, max = 3, step = 0.1)
                                ),
                                column(6,
                                       selectInput("prep_linkage", "Linkage method",
                                                   choices = c("Ward (wardu)" = "wardu",
                                                               "Average" = "average",
                                                               "Complete" = "complete",
                                                               "Single" = "single"),
                                                   selected = "wardu"),
                                       checkboxInput("prep_log2", "Log2 transform", value = TRUE),
                                       checkboxInput("prep_mean_center", "Mean‑center", value = TRUE),
                                       checkboxInput("prep_unit_variance", "Unit variance scaling", value = TRUE),
                                       checkboxInput("prep_transpose", "Transpose (pixel clustering)", value = TRUE)
                                )
                            )
                        ),
                        hr(),
                        actionButton("run_preprocess", "Run Preprocessing", class = "btn-primary", 
                                     icon = icon("play")),
                        br(), br(),
                        conditionalPanel(
                            condition = "output.preprocessing_done",
                            div(style = "text-align: center;",
                                downloadButton("download_prep_cluster", "Download cluster file (.txt)"),
                                downloadButton("download_prep_intensity", "Download intensity table (.tsv)"),
                                br(), br(),
                                actionButton("load_generated_data", "Load generated data into app", 
                                             class = "btn-success", icon = icon("upload"))
                            )
                        )
                    ),
                    box(
                        width = 8,
                        h4("Preprocessing status"),
                        verbatimTextOutput("preprocess_status"),
                        br(),
                        h4("Preview of generated data"),
                        tableOutput("prep_preview")
                    )
                )
            ),
            # ========================= SPATIAL DOMAINS =========================
            tabItem(
                tabName = "spatialdomains",
                h1("Visualization of spatial domains", align = "center"),
                box(width = 12, uiOutput("spatial_domains_doc")),
                fluidRow(
                    box(
                        width = 12,
                        column(width = 3,
                               numericInput("alpha_value_hc", "Transparency value", value = 0.5, min = 0, max = 1, step = 0.05),
                               hr(),
                               selectInput("cluster_filter", "Show only cluster:",
                                           choices = c("All" = "All"),
                                           selected = "All"),
                               hr(),
                               checkboxInput("white_bg_hc", "White background", value = FALSE)
                        ),
                        column(width = 1),
                        column(
                            width = 3,
                            textInput("label_input_hc", "Annotation label/name", placeholder = "Default"),
                            fluidRow(
                                column(width = 1),
                                actionButton("label_confirm_hc", "Confirm annotation"),
                                downloadButton("export_df_hc", "Export annotations")
                            )
                        ),
                        column(width = 3, 
                               textInput("label_input_numerichc", "Specify clusters to re-annotate separated by comma (Note: Replaces manual/lasso selection)", placeholder = NULL)
                        )
                    )
                ),
                fluidRow(
                    box(
                        width = 6,
                        plotOutput("hc_plot", width = "100%", height = "600px"),
                        downloadButton("download_hc_pdf", "Download PDF"),
                        downloadButton("download_hc_png", "Download PNG")
                    ),
                    conditionalPanel(
                        "output.imgUploaded",
                        box(width = 6, imageOutput("img_hc_plot"))
                    )
                )
            ),
            
            # ========================= INTENSITIES =========================
            tabItem(
                tabName = "intensities",
                h1("Molecule intensities", align = "center"),
                fluidRow(
                    box(
                        width = 12,
                        title = "Controls",
                        status = "primary",
                        solidHeader = TRUE,
                        column(
                            width = 3,
                            colourpicker::colourInput("mz_spot_color", "Spot Color", value = "#1f78b4"),
                            sliderInput("alpha_value_mz", "Transparency (Alpha)", min = 0, max = 1, value = 0.7, step = 0.05),
                            sliderInput("mz_contrast", "Contrast", min = 0.3, max = 5, value = 1, step = 0.05),
                            sliderInput("mz_point_size", "Point Size", min = 1, max = 10, value = 3, step = 0.5),
                            checkboxInput("white_bg_mz", "Use white background", value = FALSE),
                            actionButton("apply_mz_settings", "Apply Settings", icon = icon("refresh")),
                            br(), br(),
                            actionButton("force_gc", "Refresh", icon = icon("broom"))
                        ),
                        column(width = 1),
                        column(
                            width = 6,
                            fluidRow(
                                column(width = 6, textInput("label_input_mz", "Annotation label/name", placeholder = "Default")),
                                column(width = 6, selectizeInput("mz_selected", "Select a feature to plot", choices = NULL))
                            ),
                            fluidRow(
                                column(width = 4),
                                column(width = 4, actionButton("label_confirm_mz", "Confirm annotation")),
                                column(width = 4, downloadButton("export_df_mz", "Export annotations"))
                            )
                        )
                    )
                ),
                fluidRow(
                    box(
                        width = 12,
                        title = "Spatial Intensity Visualization",
                        status = "primary",
                        solidHeader = TRUE,
                        ggiraph::girafeOutput("mz_plot", width = "100%", height = "700px"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_mz_pdf", "Download PDF"),
                            downloadButton("download_mz_png", "Download PNG")
                        )
                    )
                ),
                conditionalPanel(
                    condition = "output.imgUploaded",
                    fluidRow(
                        box(width = 12, title = "Reference Image", imageOutput("img_mz_plot"))
                    )
                )
            ),
            
            # ========================= VIOLIN PLOT =========================
            tabItem(
                tabName = "violin",
                h1("Intensity distribution across clusters", align = "center"),
                fluidRow(
                    box(
                        width = 12,
                        column(
                            width = 4,
                            radioButtons("violin_mode", "Plot mode",
                                         choices = c("Single feature" = "single",
                                                     "All features (global distribution)" = "all"),
                                         selected = "single"),
                            radioButtons("violin_transform", "Intensity transformation",
                                         choices = c("Raw" = "raw", "Log2(x+1)" = "log2"),
                                         selected = "raw"),
                            radioButtons("violin_plot_type", "Plot type",
                                         choices = c("Violin only" = "violin",
                                                     "Boxplot only" = "boxplot",
                                                     "Violin + Boxplot" = "both"),
                                         selected = "violin"),
                            conditionalPanel(
                                condition = "input.violin_mode == 'single'",
                                selectizeInput("violin_feature", "Select metabolite", choices = NULL,
                                               options = list(placeholder = "Type to search"))
                            ),
                            conditionalPanel(
                                condition = "input.violin_mode == 'all'",
                                helpText("Shows the distribution of all intensity values (all metabolites, all pixels) within each cluster.")
                            ),
                            hr(),
                            h4("Cluster filter"),
                            selectizeInput("violin_clusters", "Filter to clusters (optional):",
                                           choices = NULL, multiple = TRUE,
                                           options = list(placeholder = "All clusters")),
                            hr(),
                            h4("Font sizes"),
                            sliderInput("violin_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
                            sliderInput("violin_axis_text_size", "Axis tick label size", min = 8, max = 20, value = 10, step = 1),
                            sliderInput("violin_title_size", "Plot title size", min = 10, max = 24, value = 14, step = 1),
                            sliderInput("violin_legend_text_size", "Legend text size", min = 8, max = 16, value = 10, step = 1),
                            hr(),
                            h4("Axis options"),
                            sliderInput("violin_x_angle", "X-axis label rotation (degrees)", 
                                        min = 0, max = 90, value = 45, step = 5),
                            checkboxInput("violin_transpose", "Transpose (clusters on Y axis)", value = FALSE)
                        ),
                        column(
                            width = 8,
                            plotlyOutput("violin_plot", height = "500px"),
                            br(),
                            div(style = "text-align: center;",
                                downloadButton("download_violin_png", "Download PNG"),
                                downloadButton("download_violin_pdf", "Download PDF")
                            )
                        )
                    )
                )
            ),
            
            # ========================= PCA PLOT =========================
            tabItem(
                tabName = "pca",
                h1("Principal Component Analysis of Metabolite Intensities", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        radioButtons("pca_transform", "Intensity transformation",
                                     choices = c("Raw" = "raw", "Log2(x+1)" = "log2"),
                                     selected = "log2"),
                        numericInput("pca_top_n", "Number of top variable metabolites", 
                                     value = 500, min = 50, max = 5000),
                        selectizeInput("pca_clusters", "Filter pixels by cluster (optional):",
                                       choices = NULL, multiple = TRUE,
                                       options = list(placeholder = "All clusters")),
                        numericInput("pca_dim_x", "PC for x‑axis", value = 1, min = 1, max = 10),
                        numericInput("pca_dim_y", "PC for y‑axis", value = 2, min = 1, max = 10),
                        sliderInput("pca_point_size", "Point size", min = 0.5, max = 3, value = 1.5, step = 0.1),
                        sliderInput("pca_point_alpha", "Point transparency", min = 0.2, max = 1, value = 0.7, step = 0.05),
                        checkboxInput("pca_add_ellipse", "Add confidence ellipses", value = TRUE),
                        hr(),
                        h4("Font sizes"),
                        sliderInput("pca_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("pca_axis_text_size", "Axis tick label size", min = 8, max = 20, value = 10, step = 1),
                        sliderInput("pca_title_size", "Plot title size", min = 10, max = 24, value = 14, step = 1),
                        sliderInput("pca_legend_text_size", "Legend text size", min = 8, max = 16, value = 10, step = 1),
                        actionButton("run_pca", "Run PCA", class = "btn-primary")
                    ),
                    box(
                        width = 9,
                        plotlyOutput("pca_plot", height = "600px"),
                        br(),
                        verbatimTextOutput("pca_variance"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_pca_png", "Download PNG"),
                            downloadButton("download_pca_pdf", "Download PDF")
                        )
                    )
                )
            ),
            
            # ========================= SCATTER PLOT =========================
            tabItem(
                tabName = "scatter",
                h1("Scatter plot of two metabolites", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        selectizeInput("scatter_x", "Metabolite for X‑axis", choices = NULL,
                                       options = list(placeholder = "Type to search")),
                        selectizeInput("scatter_y", "Metabolite for Y‑axis", choices = NULL,
                                       options = list(placeholder = "Type to search")),
                        radioButtons("scatter_transform", "Intensity transformation",
                                     choices = c("Raw" = "raw", "Log2(x+1)" = "log2"),
                                     selected = "log2"),
                        selectizeInput("scatter_clusters", "Filter pixels by cluster (optional):",
                                       choices = NULL, multiple = TRUE,
                                       options = list(placeholder = "All clusters")),
                        sliderInput("scatter_point_size", "Point size", min = 0.5, max = 5, value = 1.5, step = 0.1),
                        sliderInput("scatter_point_alpha", "Point transparency", min = 0.2, max = 1, value = 0.6, step = 0.05),
                        checkboxInput("scatter_add_line", "Add linear regression line (with confidence band)", value = FALSE),
                        hr(),
                        h4("Font sizes"),
                        sliderInput("scatter_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("scatter_axis_text_size", "Axis tick label size", min = 8, max = 20, value = 10, step = 1),
                        sliderInput("scatter_title_size", "Plot title size", min = 10, max = 24, value = 14, step = 1),
                        sliderInput("scatter_legend_text_size", "Legend text size", min = 8, max = 16, value = 10, step = 1),
                        actionButton("run_scatter", "Update plot", class = "btn-primary")
                    ),
                    box(
                        width = 9,
                        plotlyOutput("scatter_plot", height = "600px"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_scatter_png", "Download PNG"),
                            downloadButton("download_scatter_pdf", "Download PDF")
                        )
                    )
                )
            ),
            
            # ========================= UMAP PLOT =========================
            tabItem(
                tabName = "umap",
                h1("UMAP of Metabolite Intensities", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        radioButtons("umap_transform", "Intensity transformation",
                                     choices = c("Raw" = "raw", "Log2(x+1)" = "log2"),
                                     selected = "log2"),
                        numericInput("umap_top_n", "Number of top variable metabolites", 
                                     value = 1000, min = 50, max = 5000, step = 50),
                        numericInput("umap_n_neighbors", "Number of neighbors (n_neighbors)", 
                                     value = 15, min = 5, max = 100, step = 1),
                        numericInput("umap_min_dist", "Minimum distance (min_dist)", 
                                     value = 0.1, min = 0.01, max = 0.99, step = 0.05),
                        numericInput("umap_dim_x", "UMAP component for X‑axis", value = 1, min = 1, max = 5),
                        numericInput("umap_dim_y", "UMAP component for Y‑axis", value = 2, min = 1, max = 5),
                        sliderInput("umap_point_size", "Point size", min = 0.5, max = 3, value = 1.5, step = 0.1),
                        sliderInput("umap_point_alpha", "Point transparency", min = 0.2, max = 1, value = 0.7, step = 0.05),
                        checkboxInput("umap_add_ellipse", "Add confidence ellipses", value = TRUE),
                        checkboxInput("umap_swap_axes", "Swap axes (flip X and Y)", value = FALSE),
                        hr(),
                        h4("Font sizes"),
                        sliderInput("umap_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("umap_axis_text_size", "Axis tick label size", min = 8, max = 20, value = 10, step = 1),
                        sliderInput("umap_title_size", "Plot title size", min = 10, max = 24, value = 14, step = 1),
                        sliderInput("umap_legend_text_size", "Legend text size", min = 8, max = 16, value = 10, step = 1),
                        actionButton("run_umap", "Run UMAP", class = "btn-primary")
                    ),
                    box(
                        width = 9,
                        plotlyOutput("umap_plot", height = "600px"),
                        br(),
                        verbatimTextOutput("umap_info"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_umap_png", "Download PNG"),
                            downloadButton("download_umap_pdf", "Download PDF")
                        )
                    )
                )
            ),
            
            # ========================= SPATIAL STATISTICS =========================
            tabItem(
                tabName = "spatial_stats",
                h1("Spatial statistics", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        textInput("user_top_var", "Calculate statistics for this many variable features:", value = 10),
                        actionButton("test_run_spatial_stats", "Calculate")
                    ),
                    box(
                        width = 9,
                        DT::DTOutput("sp_output")
                    )
                )
            ),
            
            # ========================= SPATIAL GRADIENTS =========================
            tabItem(
                tabName = "spatial_gradients",
                h1("Spatial gradient detection", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        textInput("user_top_var_grad", "Test for spatial gradients in this many variable features:", value = 100),
                        selectInput("annotation_test_gradients", "Annotation to test",
                                    choices = c("hcdist domains" = "hc_orig",
                                                "hcdist-based annotations" = "hc_manual",
                                                "Feature-based annotations" = "mz_manual")),
                        selectizeInput("ref_group", "Reference domain", choices = NULL),
                        selectInput("summ_type", "Summarize distances with:", choices = c("Average" = "avg", "Minimum" = "min")),
                        actionButton("test_run_spatial_gradients", "Calculate")
                    ),
                    box(
                        width = 9,
                        DT::DTOutput("spgradient_output")
                    )
                )
            ),
            
            # ========================= FPCA (SPATIAL G‑FUNCTION) =========================
            tabItem(
                tabName = "fpca",
                h1("Functional PCA of Spatial G‑function Differences", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        h4("Parameters"),
                        numericInput("fpca_r_max", "Maximum distance (r max)", value = 20, step = 5),
                        numericInput("fpca_r_step", "Distance step", value = 0.5, step = 0.1),
                        numericInput("fpca_n_perms", "Number of permutations", value = 50, min = 10, step = 10),
                        numericInput("fpca_pve", "Proportion of variance explained", value = 0.99, min = 0.8, max = 1, step = 0.01),
                        numericInput("fpca_max_clusters", "Max number of clusters (by size)", value = 20, min = 3, step = 1),
                        numericInput("fpca_subsample", "Pixel subsampling (keep 1 of every N)", value = 1, min = 1, step = 1),
                        checkboxInput("fpca_use_collapsed", "Use collapsed pixels (if available)", value = TRUE),
                        hr(),
                        h4("Score Plot Appearance"),
                        sliderInput("fpca_point_size", "Point size", min = 1, max = 10, value = 3, step = 0.5),
                        sliderInput("fpca_label_size", "Cluster label size", min = 1, max = 8, value = 4, step = 0.5),
                        sliderInput("fpca_title_size", "Plot title size", min = 8, max = 24, value = 16, step = 1),
                        sliderInput("fpca_axis_title_size", "Axis title size", min = 8, max = 20, value = 14, step = 1),
                        sliderInput("fpca_axis_text_size", "Axis text size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("fpca_legend_text_size", "Legend text size", min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        h4("Variance Plot Appearance"),
                        radioButtons("fpca_line_type", "Red line represents:",
                                     choices = c("Individual variances (decreasing)" = "individual",
                                                 "Cumulative variance (increasing)" = "cumulative"),
                                     selected = "individual"),
                        sliderInput("fpca_var_title_size", "Variance plot title size", min = 8, max = 24, value = 14, step = 1),
                        sliderInput("fpca_var_axis_title_size", "Variance axis title size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("fpca_var_axis_text_size", "Variance axis text size", min = 8, max = 20, value = 10, step = 1),
                        sliderInput("fpca_var_bar_width", "Bar width", min = 0.3, max = 0.9, value = 0.7, step = 0.05),
                        hr(),
                        actionButton("run_fpca", "Run FPCA", class = "btn-primary"),
                        br(), br(),
                        downloadButton("download_fpca_scores", "Download FPCA scores")
                    ),
                    box(
                        width = 9,
                        h4("FPCA Score Plot (FPC1 vs FPC2)"),
                        plotlyOutput("fpca_score_plot", height = "500px"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_fpca_score_png", "Download PNG"),
                            downloadButton("download_fpca_score_pdf", "Download PDF")
                        ),
                        br(),
                        h4("Variance explained"),
                        plotlyOutput("fpca_var_plot", height = "300px"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_fpca_var_png", "Download PNG"),
                            downloadButton("download_fpca_var_pdf", "Download PDF")
                        ),
                        br(),
                        verbatimTextOutput("fpca_summary")
                    )
                )
            ),
            
            # ========================= DIFFERENTIAL ABUNDANCE (moved here) =========================
            tabItem(
                tabName = "diff_abund",
                h1("Differential abundance between two clusters", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        radioButtons("intensity_transform", "Intensity transformation",
                                     choices = c("Raw" = "raw", "Log2(x+1)" = "log2"),
                                     selected = "log2"),
                        selectInput("annotation_test", "Annotation column to use",
                                    choices = c("hcdist domains" = "hc_orig",
                                                "hcdist-based annotations" = "hc_manual",
                                                "Feature-based annotations" = "mz_manual"),
                                    selected = "hc_manual"),
                        selectizeInput("group1_selected", "Cluster 1 (reference)", 
                                       choices = NULL, multiple = FALSE),
                        checkboxInput("compare_to_all", 
                                      label = strong("Compare Cluster 1 vs All Other Clusters"), 
                                      value = FALSE),
                        conditionalPanel(
                            condition = "!input.compare_to_all",
                            selectizeInput("group2_selected", "Cluster 2", 
                                           choices = NULL, multiple = FALSE)
                        ),
                        selectInput("diff_test", "Test type",
                                    choices = c("Wilcoxon Rank Test" = "wilcoxon",
                                                "T-test" = "ttest",
                                                "Hellinger Distance" = "hellinger",
                                                "Spatial limma (with coordinates)" = "limma_spatial")),
                        actionButton("test_run_mz", "Run DE test", class = "btn-primary"),
                        hr(),
                        h4("Volcano plot options"),
                        conditionalPanel(
                            condition = "input.intensity_transform == 'log2'",
                            sliderInput("volcano_fc_cutoff", "log2 Fold Change cutoff",
                                        min = 0, max = 5, value = 0.1, step = 0.05)
                        ),
                        conditionalPanel(
                            condition = "input.intensity_transform == 'raw'",
                            sliderInput("volcano_fc_cutoff", "Raw fold change cutoff (difference of means)",
                                        min = 0, max = 10000, value = 100, step = 100)
                        ),
                        sliderInput("volcano_fdr_cutoff", "FDR cutoff",
                                    min = 0, max = 0.2, value = 0.05, step = 0.01),
                        selectInput("volcano_y", "Y-axis",
                                    choices = c("Test statistic" = "stat", "-log10(p-value)" = "neglogp"),
                                    selected = "stat"),
                        numericInput("volcano_top_n", "Number of top features to label",
                                     value = 20, min = 1, max = 100),
                        textInput("volcano_extra", "Extra features to label (comma-separated)", value = ""),
                        hr(),
                        h4("Font sizes"),
                        sliderInput("volcano_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("volcano_axis_text_size", "Axis tick label size", min = 8, max = 20, value = 10, step = 1),
                        sliderInput("volcano_title_size", "Plot title size", min = 10, max = 24, value = 14, step = 1),
                        sliderInput("volcano_legend_text_size", "Legend text size", min = 8, max = 16, value = 10, step = 1),
                        sliderInput("volcano_label_size", "Feature label size", min = 2, max = 6, value = 3, step = 0.2),
                        actionButton("update_volcano", "Update volcano plot", class = "btn-default")
                    ),
                    box(
                        width = 9,
                        DT::DTOutput("da_output"),
                        br(),
                        plotlyOutput("volcano_plot", height = "600px"),
                        br(),
                        div(style = "text-align: center;",
                            downloadButton("download_volcano_png", "Download PNG"),
                            downloadButton("download_volcano_pdf", "Download PDF")
                        )
                    )
                )
            ),
            
            # ========================= METABOLITE NETWORK =========================
            tabItem(
                tabName = "met_network",
                h1("Metabolite Co-localization Network", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        h4("Network Parameters"),
                        radioButtons("net_filter_mode", "Feature Selection Mode",
                                     choices = c("Global Top Variable" = "global",
                                                 "Top from Specific Domain" = "per_domain"),
                                     selected = "global"),
                        conditionalPanel(
                            "input.net_filter_mode == 'global'",
                            numericInput("net_top_n", "Number of top features (global)", 
                                         value = 300, min = 50, max = 1500, step = 50)
                        ),
                        conditionalPanel(
                            "input.net_filter_mode == 'per_domain'",
                            selectInput("net_selected_domain", "Select Domain / Cluster", choices = NULL),
                            numericInput("net_top_per_domain", "Top features from selected domain", 
                                         value = 50, min = 10, max = 200, step = 5)
                        ),
                        selectInput("net_corr_method", "Correlation method",
                                    choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                                    selected = "pearson"),
                        sliderInput("net_corr_threshold", "Absolute correlation threshold",
                                    min = 0.3, max = 0.9, value = 0.65, step = 0.01),
                        numericInput("net_min_degree", "Minimum node degree", value = 2, min = 1),
                        radioButtons("net_layout", "Layout algorithm",
                                     choices = c("Force-directed (Fruchterman-Reingold)" = "fr",
                                                 "Kamada-Kawai" = "kk",
                                                 "Circle" = "circle"),
                                     selected = "fr"),
                        checkboxInput("net_interactive", "Enable dragging", value = TRUE),
                        colourpicker::colourInput("node_color", "Node Color", value = "#1f78b4"),
                        colourpicker::colourInput("edge_color_pos", "Positive Edge Color", value = "#1f78b4"),
                        colourpicker::colourInput("edge_color_neg", "Negative Edge Color", value = "#e31a1c"),
                        
                        checkboxInput("net_show_edge_labels", "Show correlation values on edges", value = TRUE),
                        hr(),
                        h4("Plot Engine"),
                        radioButtons("net_plot_type", "Select plot type",
                                     choices = c("visNetwork" = "vis",
                                                 "plotly" = "plotly"),
                                     selected = "vis"),
                        sliderInput("net_edge_label_size", "Edge label font size", 
                                    min = 8, max = 30, value = 12, step = 1),
                        sliderInput("net_node_label_size", "Node label font size", 
                                    min = 8, max = 30, value = 14, step = 1),
                        
                        # ---- Layout saving and screenshot buttons ----
                        actionButton("save_layout", "Save Current Layout for Download", 
                                     icon = icon("save"), class = "btn-info"),
                        br(), br(),
                        
                        actionButton("run_network", "Build Network", class = "btn-primary"),
                        hr(),
                        
                        downloadButton("download_network_png", "Download PNG (static)"),
                        downloadButton("download_network_pdf", "Download PDF (static)"),
                        br(), br(),
                        
                        # NEW: screenshot button
                        actionButton("screenshot_btn", "Download Screenshot (PNG)", 
                                     icon = icon("camera"), class = "btn-success")
                    ),
                    box(
                        width = 9,
                        h4("Network Plot"),
                        uiOutput("network_plot"),        # dynamic: visNetwork or plotly
                        verbatimTextOutput("network_summary")
                    )
                )
            ),
            
            # ========================= METABOLITE OCCURRENCE =========================
            tabItem(
                tabName = "met_occurrence",
                h1("Metabolite Occurrence per Domain", align = "center"),
                fluidRow(
                    box(
                        width = 3,
                        h4("Parameters"),
                        numericInput("occ_top_n", "Number of top variable features to include",
                                     value = 100, min = 10, max = 5000, step = 10),
                        hr(),
                        h4("Statistical Test"),
                        selectInput("occ_test_method", "Test method",
                                    choices = c("Wilcoxon + FDR" = "wilcox",
                                                "Permutation (one‑sided)" = "perm"),
                                    selected = "wilcox"),
                        hr(),
                        h4("Enrichment Thresholds"),
                        sliderInput("occ_log2fc_thresh", "log2 Fold-Change threshold",
                                    min = 0.5, max = 4, value = 1, step = 0.25),
                        conditionalPanel(
                            condition = "input.occ_test_method == 'wilcox'",
                            sliderInput("occ_fdr_thresh", "FDR threshold (Benjamini‑Hochberg)",
                                        min = 0.01, max = 0.2, value = 0.05, step = 0.01)
                        ),
                        conditionalPanel(
                            condition = "input.occ_test_method == 'perm'",
                            sliderInput("occ_pval_thresh", "Permutation p‑value threshold",
                                        min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                            numericInput("occ_n_perm", "Number of permutations",
                                         value = 500, min = 100, max = 2000, step = 100)
                        ),
                        hr(),
                        h4("Plot Appearance"),
                        checkboxInput("occ_sort_clusters", "Sort bars by count (descending)", value = TRUE),
                        checkboxInput("occ_show_labels", "Show count labels on bars", value = TRUE),
                        checkboxInput("occ_show_signif", "Show significance stars (above bars)", value = TRUE),
                        sliderInput("occ_bar_alpha", "Bar transparency", min = 0.3, max = 1, value = 0.85, step = 0.05),
                        sliderInput("occ_axis_title_size", "Axis title size", min = 8, max = 20, value = 12, step = 1),
                        sliderInput("occ_axis_text_size", "Axis text size", min = 8, max = 20, value = 10, step = 1),
                        sliderInput("occ_title_size", "Plot title size", min = 10, max = 24, value = 14, step = 1),
                        sliderInput("occ_x_angle", "X-axis label rotation", min = 0, max = 90, value = 45, step = 5),
                        hr(),
                        actionButton("run_occurrence", "Calculate", class = "btn-primary"),
                        br(), br(),
                        downloadButton("download_occ_pdf", "Download PDF"),
                        downloadButton("download_occ_png", "Download PNG"),
                        hr(),
                        downloadButton("download_occ_csv", "Download table (CSV)")
                    ),
                    box(
                        width = 9,
                        plotly::plotlyOutput("occ_plot", height = "550px"),
                        br(),
                        DT::DTOutput("occ_table")
                    )
                )
            )
        )
    )
)