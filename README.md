# TBnanostring: Dimension Reduction and K-Means Clustering

An R Shiny app that applies dimension reduction (PCA, UMAP, or NMF) to the TB NanoString
gene expression data (107 genes, 179 people with active TB or latent TB infection) and runs
k-means clustering on the two reduced dimensions.

## What the app shows

- **Inputs:** dimension reduction method (PCA, UMAP, or NMF with rank 2) and number of clusters k (1–9)
- **Plot:** the two reduced dimensions, colored by k-means cluster and shaped by TB status, with cluster centers marked as X
- **Table:** cross-tabulation of clusters vs. TB status, plus the adjusted Rand index

## Files

- `app.R`: the Shiny app
- `TBnanostring.rds`: TB NanoString data (Kaipilyawar et al., *Clin Infect Dis*, 2022)

## How to run

1. Install the required R packages. NMF needs the Bioconductor package Biobase, so install that first:

   ```r
   install.packages("BiocManager")
   BiocManager::install("Biobase")
   install.packages(c("shiny", "tidyverse", "uwot", "NMF", "mclust"))
   ```

2. Download or clone this repository so that `app.R` and `TBnanostring.rds` are in the same folder.

3. Open `app.R` in RStudio and click **Run App**, or run this in R:

   ```r
   shiny::runApp("path/to/this/folder")
   ```

The app computes PCA, UMAP, and NMF once at startup, which takes a few seconds. After that,
changing the method or k updates the plot right away.