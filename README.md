
# spatialMET

A comprehensive pipeline for domain detection and annotation of spatial metabolomics data using hierarchical clustering with Shiny app integration for interactive exploration.

---

## Quick Start: Using the spatialMET Shiny App via Docker

The spatialMET Shiny application is available as a pre-built Docker container hosted on Docker Hub. To run it:

1. Ensure Docker Desktop is installed on your computer.
2. Launch the application using the command-line interface:

```bash
docker run --rm -p 3838:3838 oscareospina/spatialmet_app
```

3. Access the application by opening a web browser and navigating to:  
   [http://localhost:3838](http://localhost:3838)

---

## Building the spatialMET Shiny App from Source

To build the spatialMET Shiny application from the source code in this repository:

1. Clone and prepare the repository:
   - Download and extract this repository
   - Navigate to the repository folder using your command-line interface
2. Build the Docker container:

```bash
docker build --no-cache -t spatialmet_app .
```

3. Run the container:

```bash
docker run --rm -p 3838:3838 spatialmet_app
```

4. Access the application at:  
   [http://localhost:3838](http://localhost:3838)

---

## The spatialMET Pipeline: A Three-Step Workflow

### 1. Data Input and Preprocessing

Raw mass spectrometry imaging files are processed using [Cardinal](https://www.bioconductor.org/packages/release/bioc/html/Cardinal.html), a Bioconductor package for mass spectrometry imaging analysis. This preprocessing step converts instrument-generated files (e.g., `.ibd` and `.imzML` formats from MALDI-MSI) into structured text files containing pixel coordinates and peak intensities for each spatial location.

The resulting data formats include:

- **Intensity matrices** (peak intensities across pixels)
- **hcdist matrices** (pre-processed distance matrices for hierarchical clustering)

These formatted datasets are required inputs for the spatialMET application.

*Data preprocessing workflow converting raw MSI files to analysis-ready formats*

**Figure 1:** Overview of the data input and preprocessing workflow, showing conversion from raw instrument files to analysis-ready formats.

**File Organization:** After processing raw data and generating the required output files, transfer these processed files to the directories indicated by the red arrows in Figure 2. This organization step ensures that the downstream analysis pipeline can properly access and utilize the formatted data for spatial domain detection and visualization.

*File organization and transfer workflow showing target directories for processed data*

**Figure 2:** Workflow diagram illustrating file transfer locations. Red arrows indicate target directories where processed files should be placed to enable subsequent analysis steps.

---

### 2. Spatial Domain Detection

Spatial domains are detected using hierarchical clustering implemented through a computationally-efficient C algorithm named **hcdist**. This algorithm rapidly processes spatial metabolomics data to identify coherent tissue regions based on molecular expression patterns.

The hcdist algorithm offers:

- Computational efficiency for large MSI datasets
- Scalable performance with increasing pixel counts
- Optimal memory management for spatial clustering tasks

*Visualization of detected spatial domains using hierarchical clustering with hcdist algorithm*

**Figure 3:** Results of spatial domain detection using the hcdist hierarchical clustering algorithm, showing distinct tissue regions identified by molecular expression patterns.

---

### 3. Exploratory Data Analysis and Interactive Visualization

The spatialMET Shiny application provides comprehensive tools for downstream analysis and interactive exploration:

**Statistical Analyses:**

- Differential abundance testing to identify significant molecular differences between detected spatial domains
- Spatial autocorrelation statistics (Moran's I and Geary's C) for detecting statistically significant expression "hotspots"

**Interactive Features:**

- Manual selection and annotation of regions of interest (ROIs) for targeted investigation of specific tissue areas
- Dynamic visualization of clustering results and statistical outputs
- Export capabilities for analysis results and annotated regions

*Interactive analysis interface showing differential abundance results, spatial statistics, and ROI annotation tools*

**Figure 4:** The spatialMET Shiny application interface displaying differential abundance testing results, spatial autocorrelation statistics (Moran's I and Geary's C), and region of interest (ROI) annotation capabilities for comprehensive exploratory data analysis.

---

## Troubleshooting Common Issues

### Platform Mismatch (Apple Silicon M1/M2/M3)

If you are on an **Apple Silicon** Mac (M1, M2, M3), the pre‑compiled `hcdist` binary may be incompatible with the ARM architecture.  
To avoid “Illegal instruction” or “cannot execute binary file” errors, **build and run the image for the AMD64 platform** (x86_64 emulation):

```bash
docker build --platform linux/amd64 --no-cache -t spatialmet_app .
docker run --platform linux/amd64 --rm -p 3838:3838 spatialmet_app
```

This uses Rosetta 2 / QEMU emulation and ensures the binary runs correctly.

---

### “Permission denied” or “cannot execute binary file” errors

- Ensure the `hcdist` and `maldi_image_from_clusters` binaries are executable inside the container.  
  The Dockerfile already runs `chmod +x` on these files, but if you copied them manually, run:

```bash
chmod +x hcdist_stage1/hcdist/hcdist
chmod +x hcdist_stage1/hcdist/maldi_image_from_clusters
```

- If you see `Permission denied` for bash scripts, make sure they are executable:

```bash
chmod +x hcdist_stage1/bash_scripts/*.sh
```

---

### `dos2unix: command not found`

The script `run_pipeline.sh` may attempt to run `dos2unix` to clean line endings.  
If the utility is missing, either:

- Install it via your package manager (`brew install dos2unix` on macOS, `apt-get install dos2unix` on Linux).
- Comment out the line in `run_pipeline.sh` that calls `find ... -exec dos2unix {} +` – this is safe on Unix systems.

---

### Segmentation faults or buffer overflow

If you encounter a `*** buffer overflow detected ***` or segmentation fault, this is often due to the C code assuming a large minimum number of rows.  
Workarounds:

- **Increase buffer sizes** in the source code (see the comments in `Dockerfile` for patching).
- **Adjust filtering thresholds** in the Shiny app’s preprocessing tab (e.g., increase `Filter: present fraction` or decrease the standard deviation thresholds) so that more rows survive filtering and avoid the edge case.

---

### Out of memory (OOM) or slow performance

The clustering step can be memory‑intensive. Increase Docker’s resource limits:

- In Docker Desktop, go to **Settings** → **Resources** → increase **Memory** and **CPU**.
- You can also reduce the number of cores used in the Shiny app’s preprocessing tab.

---

### Windows path issues

If you are on Windows and using WSL2, ensure you are running Docker commands from a WSL terminal, and use forward slashes (`/`) in paths.  
If you encounter issues with line endings (CRLF), ensure your scripts are saved with Unix line endings (LF). You can convert them with:

```bash
find . -name "*.sh" -exec dos2unix {} \;
```

---

### Still stuck?

Check the logs by running the container without the `--rm` flag and inspecting the output.  
For persistent issues, please open an issue on the [GitHub repository](https://github.com/biodatalab/spatialMET) with the error message and your system details.

---

## Summary

The spatialMET pipeline offers an integrated solution for spatial metabolomics analysis, combining efficient data preprocessing with advanced clustering algorithms and interactive visualization tools. By following the three-step workflow—data preprocessing, spatial domain detection, and exploratory analysis—researchers can efficiently identify and characterize molecular patterns in tissue samples, enabling deeper insights into spatial biology and metabolomics.


