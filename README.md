# spatialMET

A pipeline for domain detection and annotation of spatial metabolomics data using hierarchical 
clustering and Shiny app integration

The _spatialMET_ pipeline consist of three steps:

1. Data input: Processing raw files with [_Cardinal_](https://www.bioconductor.org/packages/release/bioc/html/Cardinal.html).
This step raw files from equipment (e.g., _.ibd_ and _.imzML_ from MALDI-MSI) and produces
text files containing pixel coordinates and peak intensities by pixel.

2. Spatial domain detection: Detection of spatial domain using hierarchical clustering using a 
computationally-efficient C algorithm named [_hcdist_](https://gitlab.moffitt.org:8000/WelshEA/maldi-clustering-and-imaging).

3. [Exploratory data analysis](#exploratory-data-analysis): A Shiny app for user-frienly, point-and-click visualization and 
downstream analysis of hierarchical clustering results from _hcdist_. Downstream analyses 
include differential abundance test and calculation of spatial statistics for "hotspot" 
detection (Moran's I and Geary's C). The app also enables manual selection and annotation of 
regions of interest (ROIs).

## Use the [spatialMET Shiny app](#exploratory-data-analysis) Docker container hosted at Docker Hub

Make sure Docker Desktop is installed in the computer. Then, using the command-line type:

```
docker run --rm -p 3838:3838 oscareospina/spatialmet_app
```

## Building the [spatialMET Shiny app](#exploratory-data-analysis)

The spatialMET Shiny app can be build from the Dockerfile included in this repository:
1. Download this repository and de-compress it.
2. Navigate to the repository folder using the command-line
3. Build the container:
```
docker build --no-cache -t spatialmet_app .
```
4. Start the container:
```
docker run --rm -p 3838:3838 shinyapp 
```
5. Open a browser and type in the address: http://localhost:3838/

