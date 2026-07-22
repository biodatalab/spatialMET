# R/preprocess.R
# Run hcdist pipeline natively (for Linux container)

write_param_file_onecol <- function(items, filename) {
    writeLines(items, con = filename)
}
write_param_file_twocol <- function(file, cores, filename) {
    writeLines(paste(file, cores), con = filename)
}

preprocess_raw_data <- function(imzml_file, out_dir, params = list(), hcdist_root = "../hcdist_stage1") {
    defaults <- list(
        cores = max(1, parallel::detectCores() - 1),
        n_clusters = 10,                           # hard‑coded default for pixel clusters
        n_mz_clusters = 50,
        filter_present = 0.01,
        filter_unlog_sd = 250,
        filter_log2_sd = 0.5,
        minkowski_p = 1.5,
        distpow = 1.5,
        linkage = "wardu",
        log2 = TRUE,
        mean_center = TRUE,
        unit_variance = TRUE,
        transpose = TRUE
    )
    for (name in names(defaults)) {
        if (!name %in% names(params)) params[[name]] <- defaults[[name]]
    }
    
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(hcdist_root)
    
    scripts_dir <- file.path(hcdist_root, "bash_scripts")
    bin_dir <- file.path(hcdist_root, "hcdist")
    
    Sys.chmod(list.files(scripts_dir, full.names = TRUE, pattern = "\\.sh$"), "0755")
    Sys.chmod(list.files(bin_dir, full.names = TRUE, pattern = "hcdist|maldi_image_from_clusters"), "0755")
    
    link_path <- file.path(hcdist_root, "maldi-clustering-and-imaging-main")
    if (!file.exists(link_path)) {
        file.symlink(bin_dir, link_path)
    }
    
    # Set environment variables
    Sys.setenv(FILTER_PRESENT = params$filter_present)
    Sys.setenv(FILTER_UNLOG_SD = params$filter_unlog_sd)
    Sys.setenv(FILTER_LOG2_SD = params$filter_log2_sd)
    Sys.setenv(MINKOWSKI_P = params$minkowski_p)
    Sys.setenv(DISTPOW = params$distpow)
    Sys.setenv(LINKAGE = params$linkage)
    Sys.setenv(LOG2 = if (params$log2) "1" else "0")
    Sys.setenv(MEAN_CENTER = if (params$mean_center) "1" else "0")
    Sys.setenv(UNIT_VARIANCE = if (params$unit_variance) "1" else "0")
    Sys.setenv(TRANSPOSE = if (params$transpose) "1" else "0")
    
    # Step 1: Cardinal processing
    param_file1 <- file.path(out_dir, "hcdist_imzml_param_list.txt")
    write_param_file_twocol(imzml_file, params$cores, param_file1)
    system(paste("bash", file.path(scripts_dir, "go_cardinal_all.sh"), param_file1), wait = TRUE)
    
    # Locate intensity table
    base_name <- sub("\\.imzML$", "", imzml_file, ignore.case = TRUE)
    intx_file <- paste0(base_name, "_intx.tsv")
    if (!file.exists(intx_file)) {
        intx_zst <- paste0(intx_file, ".zst")
        if (file.exists(intx_zst)) {
            intx_file <- file.path(out_dir, "intensity_table.tsv")
            system(paste("zstd -d -c", intx_zst, ">", intx_file), wait = TRUE)
        } else {
            stop("Cardinal processing failed: intensity table not found.")
        }
    } else {
        file.copy(intx_file, file.path(out_dir, "intensity_table.tsv"), overwrite = TRUE)
        intx_file <- file.path(out_dir, "intensity_table.tsv")
    }
    
    # Step 2: Filter data
    intx_zst <- file.path(out_dir, "intensity_table.tsv.zst")
    system(paste("zstd -f -19", intx_file, "-o", intx_zst), wait = TRUE)
    param_file2 <- file.path(out_dir, "hcdist_filter_param_list.txt")
    write_param_file_onecol(intx_zst, param_file2)
    system(paste("bash", file.path(scripts_dir, "go_filter_data.sh"), param_file2), wait = TRUE)
    
    filtered_out <- paste0(sub("\\.zst$", "", intx_zst), "_filtered.tsv")
    if (!file.exists(filtered_out)) stop("Filtering failed.")
    filtered_final <- file.path(out_dir, "intensity_filtered.tsv")
    file.rename(filtered_out, filtered_final)
    
    # Step 3a: Feature tree
    feature_tree_out <- file.path(out_dir, "feature_tree.tre")
    system(paste("bash", file.path(scripts_dir, "go_cluster_features.sh"),
                 filtered_final, feature_tree_out, params$cores), wait = TRUE)
    if (!file.exists(feature_tree_out)) stop("Feature tree generation failed")
    
    # Step 3b: Pixel tree
    pixel_tree_out <- file.path(out_dir, "pixel_tree.tre")
    system(paste("bash", file.path(scripts_dir, "go_detect_spatial_domains.sh"),
                 filtered_final, pixel_tree_out, params$cores), wait = TRUE)
    if (!file.exists(pixel_tree_out)) stop("Pixel tree generation failed")
    
    # Step 4: Generate result files
    out_token <- file.path(out_dir, "hcdist_result")
    system(paste("bash", file.path(scripts_dir, "go_generate_hcdist_result_files.sh"),
                 filtered_final, feature_tree_out, pixel_tree_out, out_token,
                 params$n_mz_clusters, params$n_clusters), wait = TRUE)
    
    possible_files <- list.files(out_dir, pattern = paste0("^", basename(out_token), ".*"), full.names = TRUE)
    if (length(possible_files) == 0) stop("Result file generation failed.")
    
    # Find the pixel cluster file
    cluster_file <- grep("pixel.*cluster", possible_files, value = TRUE, ignore.case = TRUE)
    if (length(cluster_file) == 0) cluster_file <- grep("clusters", possible_files, value = TRUE)
    if (length(cluster_file) == 0) stop("Could not locate cluster output file.")
    
    # Use the filtered intensity file (already present) as the intensity table
    intensity_final <- filtered_final
    
    if (!file.exists(intensity_final)) stop("Filtered intensity file missing.")
    cluster_file <- cluster_file[1]
    intensity_final <- intensity_final[1]
    
    list(
        cluster_file = cluster_file,
        intensity_file = intensity_final,
        out_dir = out_dir
    )
}