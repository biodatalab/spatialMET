# spatialMET 1.0.0

## Major Updates

- **Integrated preprocessing pipeline**: The Shiny app now includes a fully functional preprocessing module that processes raw `.imzML` and `.ibd` files directly within the app. Users can now upload raw mass spectrometry imaging data and run the entire hcdist pipeline without leaving the Shiny interface.

- **Dual data loading workflow**: Two complementary pathways for loading data into the app:
  - **Upload pre-processed data** – load previously generated cluster and intensity files
  - **Preprocess raw data** – process raw `.imzML`/`.ibd` files and load results directly

- **Manual data loading control**: After preprocessing, users must explicitly click "Load generated data into app" before the data becomes available for analysis – preventing accidental automatic loading and providing greater control.

- **Sidebar visibility**: The analysis menu now appears automatically when data is loaded via either workflow (upload or preprocessing).

## New Features

- **Enhanced preprocessing parameters**: Users can now customise filtering thresholds (`filter-present`, `filter-unlog-sd`, `filter-log2-sd`), clustering parameters, and transformation options directly from the app interface.

- **Data preview**: After preprocessing, a preview of the generated cluster data is displayed, allowing users to verify successful processing before loading.

- **Download capabilities**: Processed files can be downloaded directly from the app after preprocessing.

- **Support for ZIP archives**: Users can now upload compressed `.zip` archives containing `.imzML` and `.ibd` files.

- **Standalone preprocessing script**: The `run_pipeline.sh` script is available for HPC environments and batch processing of large datasets.


## Technical Improvements

- **Docker build optimisation**: Improved Dockerfile with buffer patching and correct binary linking.
- **Symlink fixes**: Corrected symlink creation to avoid "Is a directory" errors.
- **ARM64 compatibility**: Added support for Apple Silicon (M1/M2/M3) via AMD64 emulation.
- **Standalone script documentation**: Enhanced `run_pipeline.sh` documentation for HPC and batch processing.
- **Dependency management**: Updated `dependencies.txt` with comprehensive list of required R packages.
- **Docker Hub deployment**: Pre-built image available at `yonatan2627/spatialmet_app`.

## Documentation

- Expanded README with detailed figures (Figure 1-5) illustrating the entire workflow.
- Added troubleshooting section covering common Docker issues, platform mismatches, buffer overflows, and memory management.
- Updated quick start guide with both official and custom Docker image options.
- Improved `dependencies.txt` with complete package list for local installations.
- Added comprehensive installation and usage instructions.

## Dependencies

- Added all required CRAN and Bioconductor packages to `dependencies.txt`.
- Docker image now includes all dependencies pre-installed.
- Full package list available in `dependencies.txt` for local installations.

---

## Contributors

- Yonatan Ayalew Mekonnen
- Oscar E. Ospina

---

## Links

- **Repository**: https://github.com/biodatalab/spatialMET
- **Docker Hub**: https://hub.docker.com/r/yonatan2627/spatialmet_app
- **Documentation**: See README.md for full usage instructions.
- **Issue Tracker**: https://github.com/biodatalab/spatialMET/issues

---


