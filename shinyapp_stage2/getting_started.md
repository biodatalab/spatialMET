---
title: "Getting Started"
author: "Oscar Ospina"
date: "07/09/2025"
output: 
  html_document:
    mathjax: local
    self_contained: false
---



### How to get started

1. Upload spatial domain detection results (tab-delimited format). The spatial domain results can be obtained via hierarchical clustering (e.g., [spatialMET - hcdist](https://github.com/biodatalab/spatialMET)), or other tools such as Cardinal (e.g., [spatial shrunken centroids](https://www.bioconductor.org/packages/release/bioc/vignettes/Cardinal/inst/doc/Cardinal3-stats.html)). Visit [spatialMET](https://github.com/biodatalab/spatialMET) for the required input format.
2. Upload filtered intensity matrix (metabolites in rows, pixels in columns). Visit [spatialMET](https://github.com/biodatalab/spatialMET) for the required input format.
3. Optionally, upload an accompanying tissue image for visualization purposes (PNG format). The Shiny app does not perform any modification to the image.

**AFTER UPLOADING FILES**, the visualization and analysis modules become available on the sidebar:

* **Spatial domains:** Visualize the results from spatial domain detection and, optionally, merge domains into ROIs for downstream analysis.
* **Intensities:** Visualize the pixel intensity of a metabolite. Manual ROIs can also be generated via lasso selection.
* **Differential abundance:** Use spatial domains (or manually-defined ROIs) to conduct differential abundance analysis.
* **Spatial statistics:** Calculate Moran's I and Geary's C for a series of metabolites to identify spatial patterns.
* **Spatial gradient tests:** Identify metabolites with abundance varying with respect to distance from a spatial domain or ROI.

### Instructions to create a manual annotation (ROI):
1. After loading the files, go to "Spatial domains" or "Intensities" section
2. Specify the name of the label you want to create
3. Use the blue lasso (located top-right of the plot) to select the pixels in the annotation
4. Press "Confirm annotation" to set the label
5. Repeat 2-4 until all labels are set
6. Go to "Differential abundance" to use the annotations

### Troubleshooting

* **When moving to the _Spatial domains_ tab, I receive the error: `Required columns missing in the hcdist output`. What can I do?:** 
<br/>If your data is stored in a cloud service like Dropbox, or in a network drive, please make sure that the file size is larger than zero, or make a local copy in your Desktop and try loading the data again.

