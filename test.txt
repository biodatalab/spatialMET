FROM rocker/r-ver:4.6.1

# ----------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libcurl4-openssl-dev libxml2-dev libssl-dev libuv1-dev \
    ca-certificates curl wget zstd gcc g++ make cmake pkg-config \
    libfontconfig1-dev libfreetype6-dev libcairo2-dev libharfbuzz-dev \
    libfribidi-dev libtiff5-dev libglib2.0-dev libgdal-dev \
    libudunits2-dev libproj-dev libgeos-dev libfftw3-dev libabsl-dev \
    libgsl-dev time libfreeimage-dev binutils \
    && apt-get install -y --reinstall ca-certificates openssl \
    && update-ca-certificates --fresh \
    && rm -rf /var/lib/apt/lists/*

# ----------------------------------------------------------------------
# Application directories and files
# ----------------------------------------------------------------------
RUN mkdir -p /app
COPY shinyapp_stage2 /app
COPY hcdist_stage1 /hcdist_stage1

# ----------------------------------------------------------------------
# Build hcdist and maldi_image_from_clusters from source
# ----------------------------------------------------------------------
WORKDIR /hcdist_stage1/hcdist

# Compile all common object files (including cet_colors.c for the palette)
RUN gcc -O2 -std=gnu99 -c hcdist.c text.c tree.c rand_xoshiro256.c cet_colors.c

# Strip the 'main' symbol from hcdist.o so it can be used by maldi_image_from_clusters
RUN objcopy --strip-symbol=main hcdist.o hcdist_no_main.o

# Build the main hcdist executable
RUN gcc -O2 -std=gnu99 -o hcdist hcdist.o text.o tree.o rand_xoshiro256.o cet_colors.o -lm -lpthread

# Build maldi_image_from_clusters, allowing multiple definitions
RUN gcc -O2 -std=gnu99 -o maldi_image_from_clusters maldi_image_from_clusters.c \
    hcdist_no_main.o text.o tree.o rand_xoshiro256.o cet_colors.o \
    -lm -lfreeimage -Wl,--allow-multiple-definition

# Make binaries executable
RUN chmod +x hcdist maldi_image_from_clusters

# ----------------------------------------------------------------------
# Create symlink: first remove any existing directory/symlink
# ----------------------------------------------------------------------
WORKDIR /hcdist_stage1
RUN rm -rf maldi-clustering-and-imaging-main && ln -sf /hcdist_stage1/hcdist maldi-clustering-and-imaging-main

# Return to app directory
WORKDIR /app

# Ensure all bash scripts are executable
RUN chmod +x /hcdist_stage1/bash_scripts/*.sh

# ----------------------------------------------------------------------
# R global options (repos, download method, insecure flags)
# ----------------------------------------------------------------------
RUN echo 'options(\
  repos = c(CRAN = "https://cloud.r-project.org"),\
  download.file.method = "curl",\
  timeout = 1800,\
  download.file.extra = c("-L", "--insecure", "--ssl-no-revoke")\
)' >> /usr/local/lib/R/etc/Rprofile.site

ENV BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=FALSE

# ======================================================================
# STEP 1: Bioconductor (manual repos with `devel`)
# ======================================================================

# ----------------------------------------------------------------------
# 1a. Install using devel repos (this will work with R 4.6.1)
# ----------------------------------------------------------------------
RUN Rscript -e "Sys.setenv(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS='FALSE'); \
  options(repos = c(CRAN = 'https://cloud.r-project.org', \
                    BioCsoft = 'https://bioconductor.org/packages/devel/bioc', \
                    BioCann = 'https://bioconductor.org/packages/devel/data/annotation', \
                    BioCexp = 'https://bioconductor.org/packages/devel/data/experiment', \
                    BioCworkflows = 'https://bioconductor.org/packages/devel/workflows')); \
  install.packages(c('S4Vectors','BiocGenerics','BiocParallel','DelayedArray','SummarizedExperiment','ProtGenerics','limma','Cardinal'), dependencies=TRUE, quiet=FALSE)"

# ----------------------------------------------------------------------
# 1b. Verify and re-install any missing packages (still using devel repos)
# ----------------------------------------------------------------------
RUN Rscript -e " \
  packages <- c('S4Vectors','BiocGenerics','BiocParallel','DelayedArray', \
                'SummarizedExperiment','ProtGenerics','limma','Cardinal'); \
  missing <- packages[!sapply(packages, requireNamespace, quietly=TRUE)]; \
  if (length(missing) > 0) { \
    cat('Missing Bioconductor packages:', paste(missing, collapse=', '), '\n'); \
    cat('Re-installing missing packages from devel repos...\n'); \
    options(repos = c(CRAN = 'https://cloud.r-project.org', \
                      BioCsoft = 'https://bioconductor.org/packages/devel/bioc', \
                      BioCann = 'https://bioconductor.org/packages/devel/data/annotation', \
                      BioCexp = 'https://bioconductor.org/packages/devel/data/experiment', \
                      BioCworkflows = 'https://bioconductor.org/packages/devel/workflows')); \
    install.packages(missing, dependencies=TRUE, quiet=FALSE); \
  } else { \
    cat('All Bioconductor packages already installed.\n'); \
  }"

# ----------------------------------------------------------------------
# 1c. Final verification – fail if still missing
# ----------------------------------------------------------------------
RUN Rscript -e " \
  packages <- c('S4Vectors','BiocGenerics','BiocParallel','DelayedArray', \
                'SummarizedExperiment','ProtGenerics','limma','Cardinal'); \
  missing <- packages[!sapply(packages, requireNamespace, quietly=TRUE)]; \
  if (length(missing) > 0) { \
    cat('FATAL: Still missing Bioconductor packages:', paste(missing, collapse=', '), '\n'); \
    quit(status=1); \
  } else { \
    cat('All Bioconductor packages are present.\n'); \
  }"

# ======================================================================
# STEP 2: CRAN packages (only runs if Bioconductor passes)
# ======================================================================
RUN Rscript -e "install.packages(c('shiny','httpuv','fs','cpp11','dashboardthemes','shinydashboard','shinyvalidate','DT','ggiraph','plotly','colourpicker','shinyjs','htmlwidgets','progress','sf','dplyr','ggplot2','tibble','purrr','stringr','tidyr','forcats','readr','markdown','ggnewscale','khroma','units','spdep','sfsmisc','vroom','fastcluster','future','promises','processx','data.table','magrittr','Matrix','uwot','WGCNA','irlba','visNetwork','igraph','spatstat.geom','spatstat.explore','tidyverse','ggrepel','matrixStats','mxfda','parallelly'), dependencies=TRUE, quiet=TRUE)"

# ======================================================================
# STEP 3: Explicitly install tidyverse and WGCNA (redundant for safety)
# ======================================================================
RUN Rscript -e "install.packages('tidyverse', dependencies=TRUE, quiet=FALSE)"
RUN Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org', \
                                   BioCsoft = 'https://bioconductor.org/packages/devel/bioc')); \
                install.packages(c('impute','preprocessCore','WGCNA'), \
                                 dependencies=TRUE, quiet=FALSE)"

# ======================================================================
# STEP 4: Verify critical packages (fail if missing)
# ======================================================================
RUN Rscript -e " \
  critical <- c('tidyverse','WGCNA'); \
  missing <- critical[!sapply(critical, requireNamespace, quietly=TRUE)]; \
  if (length(missing) > 0) { \
    cat('ERROR: Missing critical packages:', paste(missing, collapse=', '), '\n'); \
    quit(status=1); \
  } else { \
    cat('All critical packages are installed.\n'); \
  }"

# ----------------------------------------------------------------------
# Optional verification of all installed packages
# ----------------------------------------------------------------------
RUN Rscript -e "cat('Installed packages:\n'); print(sort(installed.packages()[,'Package']))"

# ----------------------------------------------------------------------
# Set core count and start Shiny
# ----------------------------------------------------------------------
ENV R_MAX_NUM_CORES=4
EXPOSE 3838
CMD ["R", "-e", "library(parallelly); library(parallel); shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
