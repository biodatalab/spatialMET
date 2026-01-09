##
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#

# Allow large file upload
options(shiny.maxRequestSize=10000*1024^2)

# =================================== UI =======================================

dashboardPage(
  dashboardHeader(title="spatialMET"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Import data", tabName="import", icon=icon("upload")),
      conditionalPanel("output.fileUploaded",
                       sidebarMenu(
                         menuItem("Spatial domains", tabName='spatialdomains', icon=icon('braille')),
                         menuItem("Intensities", tabName='intensities', icon=icon('bolt')),
                         menuItem("Differential abundance", tabName='diff_abund', icon=icon('chart-simple')),
                         menuItem("Spatial statistics", tabName='spatial_stats', icon=icon('map')),
                         menuItem("Spatial gradient tests", tabName='spatial_gradients', icon=icon('stairs'))
                       )
      ) # CLOSE conditionalPanel
    )
  ), # CLOSE dashboardSidebar
  dashboardBody(
    tabItems(
      tabItem(tabName="import",
              h1("spatialMET: Visualization and analysis of spatial metabolomics experiments", align="center"),
              box(width=9,
                  uiOutput("getting_started")
              ),
              box(width=3,
                  fileInput("hcdist_user_file", "Upload a hcdist result file",
                            multiple=FALSE,
                            accept=".txt"),
                  fileInput("itx_user_file", "Upload a intensity table file",
                            multiple=FALSE,
                            accept=".tsv"),
                  textInput(inputId="collapse_par",
                            label="Size of neighborhood to collapse (optional)",
                            placeholder=NULL,
                            value=NULL,
                            width='150px'),
                  fileInput("img_user_file", "Upload a tissue image (optional)",
                            multiple=FALSE,
                            accept="image/*")
              )
      ), # CLOSE import
      tabItem(tabName='spatialdomains',
              h1("Visualization of spatial domains", align="center"),
              box(width=12,
                  uiOutput("spatial_domains_doc")
              ),
              fluidRow(
                box(width=12,
                    column(width=3,
                           #sliderInput("alpha_value_hc", "Transparency value", min=0, max=1, value=0.5, step=0.1)
                           textInput(inputId="alpha_value_hc", label="Transparency value", value=0.5, placeholder=0.5)
                    ),
                    column(width=1),
                    column(width=3,
                           textInput(inputId="label_input_hc", label="Annotation label/name", placeholder="Default"),
                           fluidRow(
                             column(width=1),
                             actionButton(inputId="label_confirm_hc", label="Confirm annotation"),
                             downloadButton("export_df_hc", "Export annotations"),
                           )
                    ),
                    column(width=3,
                           textInput(inputId="label_input_numerichc", label="Specify clusters to re-annotate separated by comma (Note: Replaces manual/lasso selection)", placeholder=NULL)
                    )
                )
              ),
              fluidRow(
                box(width=6,
                    ggiraph::girafeOutput("hc_plot", width="100%", height=paste0(600, "px")),
                    downloadButton("download_hc_pdf", "Download PDF"),
                    downloadButton("download_hc_png", "Download PNG")
                ),
                conditionalPanel("output.imgUploaded",
                                 box(width=6,
                                     imageOutput("img_hc_plot")
                                 )
                )
              )
      ), # CLOSE spatialdomains
      tabItem(tabName='intensities',
              h1("Molecule intensities", align="center"),
              fluidRow(
                box(width=12,
                    column(width=4,
                           #sliderInput("alpha_value_mz", "Transparency value", min=0, max=1, value=0.5, step=0.1)
                           textInput(inputId="alpha_value_mz", label="Transparency value", value=0.5, placeholder=0.5)
                    ),
                    column(width=1),
                    column(width=6,
                           fluidRow(
                             column(width=6,
                                    textInput(inputId="label_input_mz", label="Annotation label/name", placeholder="Default")
                             ),
                             column(width=6,
                                    selectizeInput('mz_selected', "Select a feature to plot", choices=NULL)
                             )
                           ),
                           fluidRow(
                             column(width=2),
                             actionButton(inputId="label_confirm_mz", label="Confirm annotation"),
                             downloadButton("export_df_mz", "Export annotations"),
                           )
                    )
                )
              ),
              fluidRow(
                box(width=6,
                    ggiraph::girafeOutput("mz_plot", width="100%", height=paste0(600, "px")),
                    downloadButton("download_mz_pdf", "Download PDF"),
                    downloadButton("download_mz_png", "Download PNG")
                ),
                conditionalPanel("output.imgUploaded",
                                 box(width=6,
                                     imageOutput("img_mz_plot")
                                 )
                )
              )
      ), # CLOSE intensities
      tabItem(tabName='diff_abund',
              h1("Differential abundance among annotations", align="center"),
              box(width=3,
                  selectInput("annotation_test", "Annotations to test", choices=c('hcdist domains'='hc_orig',
                                                                                  'hcdist-based annotations'='hc_manual',
                                                                                  'Feature-based annotations'='mz_manual')),
                  textInput('user_itx_thr', label='Test molecules within this percentile of average intensity:', placeholder=0.1, value=0.8),
                  selectizeInput('group1_selected', "Cluster 1", choices=NULL),
                  selectizeInput('group2_selected', "Cluster 2", choices=NULL),
                  selectInput("diff_test", "Test type", choices=c('Wilcoxon Rank Test'='wilcoxon', 'Hellinger Distance'='hellinger', 'T-test'='ttest')),
                  fluidRow(
                    column(width=2,
                           actionButton(inputId="test_run_mz", label="Run tests")
                    )
                  ),
                  br(),
                  fluidRow(
                    column(width=2)
                  ),
              ),
              box(width=9,
                  column(width=12,
                         DT::DTOutput("da_output")
                  )
              )
      ), # CLOSE diff_abund
      tabItem(tabName='spatial_stats',
              h1("Spatial statistics", align="center"),
              box(width=3,
                  textInput('user_top_var', label='Calculate statistics for this many variable features:', placeholder=10, value=10),
                  fluidRow(
                    column(width=2,
                           actionButton(inputId="test_run_spatial_stats", label="Calculate")
                    )
                  ),
                  br(),
                  fluidRow(
                    column(width=2)
                  ),
              ),
              box(width=9,
                  column(width=12,
                         DT::DTOutput("sp_output")
                  )
              )
      ), # CLOSE spatial_stats
      tabItem(tabName='spatial_gradients',
              h1("Spatial gradient detection", align="center"),
              box(width=3,
                  textInput('user_top_var_grad', label='Test for spatial gradients in this many variable features:', placeholder=10, value=10),
                  selectInput("annotation_test_gradients", "Annotation to test", choices=c('hcdist domains'='hc_orig',
                                                                                           'hcdist-based annotations'='hc_manual',
                                                                                           'Feature-based annotations'='mz_manual')),
                  selectizeInput('ref_group', "Reference domain", choices=NULL),
                  selectInput("summ_type", "Summarize distances with:", choices=c('Average'='avg', 'Minimum'='min')),
                  fluidRow(
                    column(width=2,
                           actionButton(inputId="test_run_spatial_gradients", label="Calculate")
                    )
                  ),
                  br(),
                  fluidRow(
                    column(width=2)
                  ),
              ),
              box(width=9,
                  column(width=12,
                         DT::DTOutput("spgradient_output")
                  )
              )
      ) # CLOSE spatial_gradients
    ) # CLOSE ITEMS
  ) # CLOSE DASHBOARD
)

