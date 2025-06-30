options(shiny.maxRequestSize = 720*1024^2)
library(shiny)
library(shinyjs)
library(fgsea)
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(plotly)
library(openxlsx)
library(ComplexHeatmap)
library(pheatmap)
library(ggrepel)

# Base data folder
data_dir <- "./data"
global_example_path <- file.path(data_dir, "global", "example.csv")
global_mapping_path <- file.path(data_dir, "global", "mapping.csv")
global_stitch_path  <- file.path(data_dir, "global", "stitch.tsv")


# Find available species folders automatically
species_folders <- list.dirs(file.path(data_dir, "species"), full.names = FALSE, recursive = FALSE)

# UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Metabolite Pathway Enrichment Explorer"),
  
  # Toggle Sidebar Button
  actionButton("toggleSidebar", "Toggle Sidebar", icon = icon("arrows-alt-h")),
  tags$hr(),
  
  fluidRow(
    column(
      width = 3,
      id = "sidebarPanel",  # Give an ID to toggle this
      wellPanel(
        selectInput("preset_species", "Choose Preset Species (or Custom Upload):",
                    choices = c("Custom Upload", species_folders),
                    selected = "Custom Upload"),
        conditionalPanel(
          condition = "input.preset_species == 'Custom Upload'",
          fileInput("pathway_file", "Upload Pathway vs Metabolites CSV", accept = ".csv"),
          fileInput("mapping_file", "Upload Mapping File (Optional)", accept = ".csv"),
          fileInput("stitch_file", "Upload STITCH Interaction File (Optional)", accept = c(".tsv", ".txt"))
        ),
        fileInput("example_file", "Upload Example Data CSV", accept = ".csv"),
        textAreaInput("metabolites", "Paste Input Metabolite IDs (one per line):", value = "", rows = 5),
        numericInput("top_n", "Top N Pathways to Show", value = 10, min = 1),
        numericInput("pval_cutoff", "Adjusted p-value Cutoff", value = 1, step = 0.01),
        actionButton("run", "Run Enrichment Analysis")
      )
    ),
    
    column(
      width = 9,
      tabsetPanel(
        tabPanel("Pathway Plot", plotlyOutput("pathwayPlot")),
        tabPanel("Impact Plot", plotlyOutput("impactPlot")),
        tabPanel("Metabolite Set Enrichment", plotlyOutput("gseaPlot")),
        tabPanel("Metabolite Centrality", plotlyOutput("rbcPlot")),
        tabPanel("Membership Plot", plotOutput("membershipPlot")),
        tabPanel("Network Plot", plotOutput("networkPlot")),
        tabPanel("Interaction Plot", plotOutput("interactionPlot")),
        tabPanel("Heatmap", plotOutput("heatmapPlot"))
      )
    )
  )
)

enrichmet <- function(inputMetabolites, PathwayVsMetabolites, example_data, top_n = 100, p_value_cutoff = 1, kegg_lookup = NULL, mapping_df = NULL, stitch_df = NULL) {
 
  matrix_to_list <- function(pws) {
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }
  
  prepare_gmt <- function(gmt_file, metabolites_in_data, savefile = FALSE) {
    gmt <- gmtPathways(gmt_file)
    hidden <- unique(unlist(gmt))
    mat <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt))
    for (i in 1:dim(mat)[2]) {
      mat[, i] <- as.numeric(hidden %in% gmt[[i]])
    }
    hidden1 <- intersect(metabolites_in_data, hidden)
    mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1, ]) > 5)]]
    final_list <- matrix_to_list(mat)
    if (savefile) {
      saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
    }
    return(final_list)
  }
  
  PathwayVsMetabolites$description <- "https://www.genome.jp/kegg/pathway.html#metabolism"
  
  convert_to_gmt <- function(pathway, description, metabolites) {
    pathway_underscore <- gsub(" ", "_", pathway)
    metabolites_vector <- unlist(strsplit(metabolites, ","))
    gmt_line <- paste(pathway_underscore, description, paste(metabolites_vector, collapse = "\t"), sep = "\t")
    return(gmt_line)
  }
  
  gmt_data <- mapply(convert_to_gmt, PathwayVsMetabolites$Pathway, PathwayVsMetabolites$description, PathwayVsMetabolites$Metabolites, SIMPLIFY = TRUE)
  gmt_file <- "output.gmt"
  writeLines(gmt_data, gmt_file)
  
  data <- PathwayVsMetabolites %>%
    mutate(Metabolites = strsplit(Metabolites, ",")) %>%
    unnest(Metabolites)
  
  allMetabolitesSet <- unique(data$Metabolites)
  
  # Centrality
  edge_list_metabolites <- data.frame(from = unlist(data$Metabolites), to = rep(data$Pathway, lengths(data$Metabolites)))
  g_metabolites <- graph_from_data_frame(d = edge_list_metabolites, directed = FALSE)
  betweenness_metabolites <- betweenness(g_metabolites, directed = FALSE, normalized = TRUE)
  metabolite_centrality <- data.frame(Metabolite = names(betweenness_metabolites), RBC_Metabolite = betweenness_metabolites)
  
  input_metabolite_centrality <- metabolite_centrality %>%
    filter(Metabolite %in% inputMetabolites) %>%
    arrange(desc(RBC_Metabolite))
  
  # Pathway enrichment
  results <- list()
  for (i in 1:nrow(PathwayVsMetabolites)) {
    row <- PathwayVsMetabolites[i, ]
    pathway <- row$Pathway
    pathwayMetabolites <- unlist(strsplit(row$Metabolites, ","))
    matchedMet <- intersect(pathwayMetabolites, inputMetabolites)
    
    a <- length(matchedMet)
    b <- length(setdiff(inputMetabolites, pathwayMetabolites))
    c <- length(setdiff(pathwayMetabolites, inputMetabolites))
    d <- length(setdiff(allMetabolitesSet, union(inputMetabolites, pathwayMetabolites)))
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    fisher_test_result <- fisher.test(contingency_table, alternative = "two.sided")
    
    matched_centrality <- metabolite_centrality %>%
      filter(Metabolite %in% matchedMet)
    
    all_centrality <- metabolite_centrality %>%
      filter(Metabolite %in% pathwayMetabolites)
    
    impact <- ifelse(nrow(all_centrality) > 0, sum(matched_centrality$RBC_Metabolite) / sum(all_centrality$RBC_Metabolite), 0)
    coverage <- length(matchedMet) / length(pathwayMetabolites)
    
    results[[i]] <- list(
      Pathway = pathway,
      P_value = fisher_test_result$p.value,
      Log_P_value = -log10(fisher_test_result$p.value),
      Impact = impact,
      Coverage = coverage
    )
  }
  
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  results_df$Adjusted_P_value <- p.adjust(results_df$P_value, method = "BH")
  
  significant_results_df <- results_df %>%
    filter(Adjusted_P_value < p_value_cutoff) %>%
    arrange(desc(Log_P_value))
  
  if (!is.null(top_n)) {
    significant_results_df <- head(significant_results_df, top_n)
  }
  
 
  if (nrow(significant_results_df) > 0) {
    # Preprocess: order Pathway levels outside ggplot
    significant_results_df <- significant_results_df %>%
      mutate(
        PathwayLabel = factor(Pathway, levels = Pathway[order(Log_P_value, decreasing = FALSE)]),
        Tooltip = paste0("Pathway: ", Pathway, "\n",
                         "-log10(P): ", round(Log_P_value, 2), "\n",
                         "Adj P-value: ", signif(Adjusted_P_value, 3))
      )
    
    
    # Updated Pathway Plot: Bubble plot sorted by -logPvalue
    pathway_plot <- ggplot(significant_results_df,
                           aes(x = PathwayLabel,
                               y = Log_P_value,
                               size = Log_P_value,
                               fill = Adjusted_P_value,
                               text = Tooltip)) +  # custom tooltip
      geom_point(shape = 21, color = NA, alpha = 0.7) +
      scale_size_continuous(range = c(3, 10), guide = "none") +
      scale_fill_gradient(
        low = "red", high = "blue",
        breaks = c(0.05, 0.5, 1),
        labels = c("0", "0.05", "1")
      ) +
      labs(x = "Pathway", y = "-log10(P-value)", fill = "Adj P-value") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      )

    # tooltip column
    significant_results_df <- significant_results_df %>%
      mutate(
        ImpactTooltip = paste0("Pathway: ", Pathway, "\n",
                               "Impact: ", round(Impact, 3), "\n",
                               "-log10(P): ", round(Log_P_value, 2), "\n",
                               "Raw P-value: ", signif(P_value, 3))
      )
    
    impact_plot <- ggplot(significant_results_df,
                          aes(x = Impact,
                              y = Log_P_value,
                              size = Impact,
                              color = P_value,
                              text = ImpactTooltip)) +  # For plotly
      geom_point(alpha = 0.85) +
      geom_text_repel(aes(label = Pathway), size = 3, max.overlaps = 100,
                      box.padding = 0.3, show.legend = FALSE) +
      scale_size_continuous(range = c(3, 10)) +
      scale_color_gradient(low = "red", high = "blue", name = "P-value") +
      labs(
        title = "Pathway Impact vs Significance",
        x = "Pathway Impact",
        y = "-log10(P-value)",
        size = "Impact"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
  } else {
    warning("No pathways passed the p-value cutoff.")
    pathway_plot <- NULL
    impact_plot <- NULL
  }
 
  # GSEA Plot
  example_filtered_data <- example_data %>%
    arrange(pval) %>%
    distinct(met_id, .keep_all = TRUE) %>%
    filter(met_id != "No Metabolites found")
  
  meta <- example_filtered_data$met_id
  bg_metabolites <- prepare_gmt(gmt_file, meta, savefile = FALSE)
  rankings <- sign(example_filtered_data$log2fc) * (-log10(as.numeric(example_filtered_data$pval)))
  names(rankings) <- example_filtered_data$met_id
  rankings <- sort(rankings, decreasing = TRUE)
  
  MSEAres <- fgsea(pathways = bg_metabolites, stats = rankings, scoreType = 'std', minSize = 10, maxSize = 500, nproc = 1)
  MSEAres$input_count <- sapply(MSEAres$leadingEdge, length)
  
  # Sort the GSEA results by p-value (ascending order)
  MSEAres <- MSEAres %>%
    arrange(pval)  # Sort by p-value in ascending order
  
  # Reversing the factor levels to make sure the lowest p-values (most significant) are at the top
  MSEAres$pathway <- factor(MSEAres$pathway, levels = rev(MSEAres$pathway))
  
  # Now create the gsea_plot
  gsea_plot <- ggplot(MSEAres, aes(x = -log10(pval), y = pathway, size = input_count, color = NES)) +
    geom_point() +
    labs(x = "-log10(p-value)", y = "Pathway", size = "Metabolite count", color = "NES") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal()
  
  write.xlsx(significant_results_df, "pathway_enrichment_results.xlsx")
  write.xlsx(MSEAres, "gsea_results.xlsx")
  
  # Centrality Plot (RBC Plot)
  if (!is.null(kegg_lookup)) {
    input_metabolite_centrality <- input_metabolite_centrality %>%
      left_join(kegg_lookup, by = c("Metabolite" = "kegg_id")) %>%
      mutate(Metabolite = coalesce(name, Metabolite))
  }
  
  rbc_plot <- ggplot(input_metabolite_centrality, aes(x = reorder(Metabolite, RBC_Metabolite), y = RBC_Metabolite, fill = RBC_Metabolite)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "orange", high = "red3") +
    labs(x = "Metabolite", y = "Relative Betweenness Centrality", fill = "RBC") +
    theme_minimal()
  
  # Network Graph (Metabolite-Pathway Network)
  df <- PathwayVsMetabolites %>%
    mutate(Metabolite = strsplit(Metabolites, ",")) %>%
    unnest(Metabolite) %>%
    filter(Metabolite %in% inputMetabolites)
  
  edge_list_all <- PathwayVsMetabolites %>%
    mutate(Metabolite = strsplit(Metabolites, ",")) %>%
    unnest(Metabolite) %>%
    select(Pathway, Metabolite)
  
  g_all <- graph_from_data_frame(edge_list_all, directed = FALSE)
  bet_all <- betweenness(g_all, directed = FALSE, normalized = TRUE)
  bet_df <- tibble(name = names(bet_all), centrality = as.numeric(bet_all))
  
  edges <- df %>% select(Pathway, Metabolite)
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  V(g)$type <- ifelse(V(g)$name %in% df$Pathway, "Pathway", "Metabolite")
  V(g)$centrality <- NA
  V(g)$centrality[V(g)$type == "Metabolite"] <- bet_df$centrality[match(V(g)$name[V(g)$type == "Metabolite"], bet_df$name)]
  V(g)$group <- ifelse(V(g)$type == "Pathway", V(g)$name, df$Pathway[match(V(g)$name, df$Metabolite)])
  
  if (!is.null(kegg_lookup)) {
    V(g)$name <- ifelse(V(g)$type == "Metabolite",
                        kegg_lookup$name[match(V(g)$name, kegg_lookup$kegg_id)],
                        V(g)$name)
  }
  
  network_layout <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +
    geom_node_point(aes(color = centrality, size = centrality)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(3, 8)) +
    theme_minimal() +
    labs(title = "Metabolite-Pathway Network", color = "Centrality")
  
  network_plot <- network_layout + theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
  
  # Enrichment Heatmap
  enriched_pathways <- significant_results_df$Pathway
  data_filtered <- data %>%
    filter(Pathway %in% enriched_pathways, Metabolites %in% inputMetabolites)
  
  # Convert KEGG IDs to names for display if kegg_lookup is provided
  if (!is.null(kegg_lookup)) {
    data_filtered <- data_filtered %>%
      left_join(kegg_lookup, by = c("Metabolites" = "kegg_id")) %>%
      mutate(Metabolites = coalesce(name, Metabolites))
  }
  
  heatmap_matrix <- table(data_filtered$Metabolites, data_filtered$Pathway)
  heatmap_matrix <- as.matrix(heatmap_matrix)
  
  logp_vec <- significant_results_df$Log_P_value
  names(logp_vec) <- significant_results_df$Pathway
  
  heatmap_values <- sweep(heatmap_matrix, 2, logp_vec[colnames(heatmap_matrix)], `*`)
  # Render heatmap with modified legend
  # Calculate appropriate breaks
  heatmap_plot <- Heatmap(
    heatmap_values,
    name = "-logP-value",
    col = circlize::colorRamp2(
      c(0, 3, 6),
      c("white", "blue", "red")
    ),
    cluster_columns = TRUE,
    border=TRUE,
    column_names_gp = gpar(fontsize = 9),  # X-axis label size
    row_names_gp = gpar(fontsize = 9)
  )
  # === Pathway Membership Plot ===
  if (!is.null(PathwayVsMetabolites)) {
    long_pathway_df <- PathwayVsMetabolites %>%
      tidyr::separate_rows(Metabolites, sep = ",") %>%
      dplyr::rename(pathway_name = Pathway, metabolite = Metabolites)
    
    matched_df <- long_pathway_df %>%
      dplyr::filter(metabolite %in% inputMetabolites)
    
    if (nrow(matched_df) > 0) {
      # Convert KEGG IDs to names if kegg_lookup is provided
      if (!is.null(kegg_lookup)) {
        matched_df <- matched_df %>%
          dplyr::left_join(kegg_lookup, by = c("metabolite" = "kegg_id")) %>%
          dplyr::mutate(name = ifelse(is.na(name), metabolite, name))
      } else {
        matched_df <- matched_df %>% dplyr::mutate(name = metabolite)
      }
      
      membership_matrix <- table(matched_df$name, matched_df$pathway_name)
      membership_matrix <- as.matrix(membership_matrix)
      membership_matrix[membership_matrix > 1] <- 1
      
      n_rows <- nrow(membership_matrix)
      n_cols <- ncol(membership_matrix)
      
      # Dynamically set size based on number of rows/cols
      heatmap_width <- unit(0.3 * n_cols, "cm")
      heatmap_height <- unit(0.4 * n_rows, "cm")
      
      membership_plot <- ComplexHeatmap::Heatmap(
        membership_matrix,
        name = "Membership",
        col = c("0" = "white", "1" = "steelblue"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_title = "Metabolites",
        column_title = "Pathways",
        heatmap_legend_param = list(title = "Member"),
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 8),
        width = heatmap_width,
        height = heatmap_height,
        border = TRUE
      )
    } else {
      warning("No matching input metabolites found in PathwayVsMetabolites.")
      membership_plot <- NULL
    }
  } else {
    warning("PathwayVsMetabolites not provided.")
    membership_plot <- NULL
  }
  # Interaction plot
  interaction_plot <- NULL  # Default if no data provided
  
  # Optional interaction plot if mapping_df and stitch_df are provided
  if (!is.null(mapping_df) && !is.null(stitch_df)) {
    
    # Filter mapping_df to only inputMetabolites (KEGG IDs)
    vertex_df <- mapping_df %>%
      filter(KEGG_ID %in% inputMetabolites) %>%
      filter(!is.na(PubChem_CID)) %>%
      distinct(KEGG_ID, .keep_all = TRUE)
    
    # Merge with KEGG lookup (if available) to get names
    if (!is.null(kegg_lookup)) {
      vertex_df <- vertex_df %>%
        left_join(kegg_lookup, by = c("KEGG_ID" = "kegg_id")) %>%
        mutate(display_name = coalesce(name, KEGG_ID))
    } else {
      vertex_df <- vertex_df %>%
        mutate(display_name = KEGG_ID)
    }
    
    vertex_df <- vertex_df %>%
      mutate(display_name = stringr::str_trunc(display_name, 25)) %>%
      select(STITCH_ID, everything())
    
    # Filter STITCH edges to valid vertices
    valid_edges <- stitch_df %>%
      filter(combined_score >= 100) %>%
      semi_join(vertex_df, by = c("chemical1" = "STITCH_ID")) %>%
      semi_join(vertex_df, by = c("chemical2" = "STITCH_ID")) %>%
      distinct(chemical1, chemical2, .keep_all = TRUE)
    
    # Build graph
    g <- igraph::graph_from_data_frame(
      d = valid_edges %>% mutate(weight = scales::rescale(combined_score, to = c(0.1, 1))),
      directed = FALSE,
      vertices = vertex_df %>%
        filter(STITCH_ID %in% c(valid_edges$chemical1, valid_edges$chemical2))
    )
    
    # Graph metrics
    igraph::V(g)$degree <- igraph::degree(g)
    igraph::V(g)$betweenness <- igraph::betweenness(g)
    igraph::V(g)$component <- igraph::components(g)$membership
    
    # Plot
    set.seed(42)
    interaction_plot <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(aes(alpha = weight), color = "grey50") +
      ggraph::geom_node_point(aes(size = degree, color = as.factor(component)), alpha = 0.8) +
      ggraph::geom_node_text(aes(label = stringr::str_wrap(display_name, 12)),
                             size = 3, repel = TRUE, max.overlaps = 20) +
      scale_color_discrete(name = "Network Component") +
      scale_size_continuous(name = "Degree", range = c(3, 10)) +
      labs(title = "Metabolite Interaction Network",
           subtitle = paste(igraph::vcount(g), "compounds with", igraph::ecount(g), "interactions")) +
      ggraph::theme_graph() +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  
  # Return everything
  return(list(
    pathway_plot = pathway_plot,
    impact_plot = impact_plot,
    gsea_plot = gsea_plot,
    rbc_plot = rbc_plot,
    network_plot = network_plot,
    heatmap_plot = heatmap_plot,
    membership_plot = membership_plot,
    interaction_plot=interaction_plot
  ))
}

server <- function(input, output,session) {
  
  observeEvent(input$toggleSidebar, {
    runjs("
    $('#sidebarPanel').toggle();
    if ($('#sidebarPanel').is(':visible')) {
      $('.col-sm-9').removeClass('col-sm-12').addClass('col-sm-9');
    } else {
      $('.col-sm-9').removeClass('col-sm-9').addClass('col-sm-12');
    }
  ")
  })
  
  
  # Preload global data once
  
  global_mapping_df <- read.csv(global_mapping_path)
  global_stitch_df <- read.table(global_stitch_path, header = TRUE, sep = "\t")
  
  results <- eventReactive(input$run, {
    # Determine which pathway/metabolite dataset to load
    if (input$preset_species != "Custom Upload") {
      # Preloaded species
      req(input$example_file)
      pathway_file_path <- file.path(data_dir, "species", input$preset_species, "pathway.csv")
      PathwayVsMetabolites <- read.csv(pathway_file_path, stringsAsFactors = FALSE)
      example_data <- read.csv(input$example_file$datapath)
      mapping_df <- global_mapping_df
      stitch_df <- global_stitch_df
    } else {
      req(input$pathway_file, input$example_file)
      PathwayVsMetabolites <- read.csv(input$pathway_file$datapath, stringsAsFactors = FALSE)
      example_data <- read.csv(input$example_file$datapath)
      # Handle optional files:
      mapping_df <- if (!is.null(input$mapping_file)) {
        read.csv(input$mapping_file$datapath)
      } else {
        NULL
      }
      
      stitch_df <- if (!is.null(input$stitch_file)) {
        read_tsv(input$stitch_file$datapath)
      } else {
        NULL
      }
    }
  
    inputMetabolites <- strsplit(input$metabolites, "\n")[[1]] %>% trimws() %>% na.omit()
    
    
    enrichmet(inputMetabolites, PathwayVsMetabolites, example_data,mapping_df = mapping_df,stitch_df = stitch_df,
              top_n = input$top_n, p_value_cutoff = input$pval_cutoff)
  })
  
  output$pathwayPlot <- renderPlotly({ req(results()); ggplotly(results()$pathway_plot,tooltip = "text")})
  output$impactPlot  <- renderPlotly({ req(results()) ;ggplotly(results()$impact_plot,tooltip = "text") })
  output$gseaPlot    <- renderPlotly({ req(results()) ;ggplotly(results()$gsea_plot) })
  output$rbcPlot     <- renderPlotly({ req(results()) ;ggplotly(results()$rbc_plot) })
  output$membershipPlot <- renderPlot({ req(results()) ;results()$membership_plot },height = 800)
  output$networkPlot <- renderPlot({ req(results()) ;results()$network_plot },height = 800)
  
  output$interactionPlot <- renderPlot({
    req(results())
    if (!is.null(results()$interaction_plot)) {
      results()$interaction_plot
    }
  })
  
  output$heatmapPlot <- renderPlot({
    req(results())
    if (!is.null(results()$heatmap_plot)) {
      grid::grid.newpage()
      ComplexHeatmap::draw(results()$heatmap_plot)
    }
  },height = 700)
}

shinyApp(ui = ui, server = server)
