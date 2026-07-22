---
title: "Getting Started with spatialMET Shiny App"
author: "Yonatan Ayalew Mekonnen"
date: "06/19/2026"
output: 
  html_document:
    mathjax: local
    self_contained: false
---

### Overview

Welcome to **spatialMET** – an interactive Shiny application for **visualization and analysis of spatial metabolomics** data.  
This tool builds on the [spatialMET R package](https://github.com/biodatalab/spatialMET) and provides a user‑friendly interface to:

- Explore spatial domain maps and metabolite intensities
- Perform differential abundance testing between domains
- Run spatial statistics (Moran's I, Geary's C) and gradient tests
- Apply dimensionality reduction (PCA, UMAP) and visualise clusters
- Generate co‑localisation networks and domain‑enriched metabolite profiles
- Use functional PCA (FPCA) on spatial G‑functions to compare domain spatial patterns

All analyses are **parallelised** to handle large datasets efficiently.

---

### How to get started

1. **Upload spatial domain detection results** (tab‑delimited `.txt` file).  
   - This file must contain at least the columns: `Pixel`, `X`, `Y`, `cluster`, `color` (with cluster labels).  
   - Optional columns: `tree order`, `cluster_size`.  
   - If the file has only 4 columns (Pixel, X, Y, cluster), the app will automatically assign colours.  
   - The expected format matches the output of `spatialMET::hcdist()` or the `spatial shrunken centroids` from `Cardinal`.  
   - **Example** (first few lines): 

     ```
     Pixel	X	Y	tree order	cluster	cluster_size	color
     pixel_1	10	20	1	1	150	#FF0000
     pixel_2	12	22	2	1	150	#FF0000
     ```

- See the [spatialMET documentation](https://github.com/biodatalab/spatialMET) for details.

2. **Upload the filtered intensity matrix** (tab‑delimited `.tsv`).  
- **Rows** = metabolites/features (e.g., m/z values).  
- **Columns** = pixels (must match the `Pixel` IDs in the domain file).  
- The first column should contain feature names, and subsequent columns are pixel intensities (numeric).  
- **Example**:

     ```
     mz_val        pixel_1  pixel_2  pixel_3  ...
     mz_153.0187   0.25     0.00     0.12     ...
     mz_211.0007   1.20     0.80     1.50     ...
     ```

- Large matrices (thousands of metabolites, tens of thousands of pixels) are handled using `DelayedArray` and parallel processing.

3. *(Optional)* **Upload a tissue image** (PNG format) – this will be displayed alongside the spatial plots for reference. The app does **not** modify or align the image; it is for visual guidance only.

> **After uploading all three files**, the sidebar menu will **unlock** all analysis modules.

---

### Modules in detail

| Module | Description |
|--------|-------------|
| **Spatial domains** | Visualise cluster assignments with optional transparency, background toggle, and cluster‑specific filtering. You can manually re‑annotate clusters (merge or rename) via the lasso tool. |
| **Intensities** | Plot the spatial distribution of a selected metabolite. Adjust contrast, point size, and colour. You can also create manual ROIs (Regions of Interest) using the lasso and give them a label. |
| **Violin plot** | Compare intensity distributions across clusters (single feature or all features combined). Supports log‑transformation, boxplot/violin overlay, and cluster filtering. |
| **PCA** | Principal Component Analysis on the most variable metabolites. Filter by cluster, choose PC axes, and add confidence ellipses. Uses fast `irlba` for speed. |
| **Scatter plot** | Correlation between two metabolites, coloured by cluster. Add a linear regression line. |
| **UMAP** | Non‑linear dimensionality reduction with configurable neighbours and minimum distance. Also supports cluster colouring and ellipses. |
| **Spatial statistics** | Compute **Moran's I** and **Geary's C** for the top N variable metabolites to quantify spatial autocorrelation. |
| **Spatial gradient tests** | Identify metabolites whose abundance changes with distance from a reference domain. Uses linear models and Spearman correlation. |
| **FPCA (Spatial G‑function)** | Functional PCA on the difference between observed and permuted G‑functions (nearest‑neighbour distance distributions). Compares spatial patterns across clusters. |
| **Differential abundance** | Compare metabolite intensities between two clusters (or one cluster vs all others) using Wilcoxon, t‑test, Hellinger distance, or a spatial limma model that accounts for pixel coordinates. Results include a volcano plot and table. |
| **Metabolite Network** | Build a correlation network from the most variable features (either globally or within a selected domain). Interactive network visualisation with adjustable thresholds. |
| **Metabolite Occurrence** | Identify metabolites that are preferentially enriched in each domain using either Wilcoxon + FDR or a permutation‑based test. Results are shown as bar plots and a sortable table. |

---

### Creating manual annotations (ROIs)

1. Go to **Spatial domains** or **Intensities** tab.
2. In the **"Annotation label/name"** text box, type a name for your ROI (e.g., "Tumour_region").
3. Use the **blue lasso tool** (top‑right of the plot) to select pixels belonging to that region.
4. Click **"Confirm annotation"** – the selected pixels will be assigned that label.
5. Repeat steps 2‑4 for other ROIs.
6. Once all ROIs are defined, switch to the **Differential abundance** tab – your custom annotations will appear as options for group selection.

> **Note:** Lasso selection is available only in the **interactive** (Plotly) plots. For static plots, use the download buttons.

---

### Performance tips

- For large datasets (> 50,000 pixels), the app will **automatically subsample** pixels (keep 1 out of every 10) to keep the interface responsive. You can adjust the **"Pixel subsampling"** parameter in the FPCA section.
- The **collapse pixels** option (via the `collapse_par` input) merges adjacent pixels into blocks, reducing the data size dramatically and speeding up downstream analyses.
- Most heavy calculations are **parallelised** using all available CPU cores. If you experience memory issues, reduce the number of features (e.g., in PCA, UMAP, or network modules).
- Use the **"Refresh"** button (in the Intensities tab) to free up memory after large computations.

---

### Troubleshooting

| Problem | Solution |
|---------|----------|
| **Error: `Required columns missing in the hcdist output`** | Ensure your domain file has at least columns: `Pixel`, `X`, `Y`, `cluster`. If the file is empty or stored on a cloud drive, download a local copy and try again. |
| **The app freezes or takes too long** | Reduce the number of features (e.g., in PCA/UMAP set a smaller top‑N). Use the **collapse pixels** option to reduce pixel count. Close other browser tabs to free memory. |
| **I see a warning about `is.na()` applied to expression** | This is a harmless warning from `plotly` when converting ggplot labels with mathematical expressions. Ignore it; the plot still works. |
| **Volcano plot labels overlap** | Adjust the **"Number of top features to label"** or manually specify features in the "Extra features to label" field. The `ggrepel` algorithm will try to avoid overlaps in the static version (download PNG/PDF). |
| **Memory errors with spatial statistics** | The app uses row‑wise parallel processing to minimise memory. If you still get errors, reduce the number of features or increase the `future.globals.maxSize` option (see R console for instructions). |

---

### Further information

- For a detailed tutorial and examples, visit the [spatialMET GitHub repository](https://github.com/biodatalab/spatialMET).
- For questions or bug reports, please open an issue on GitHub.

**Enjoy exploring your spatial metabolomics data!**



