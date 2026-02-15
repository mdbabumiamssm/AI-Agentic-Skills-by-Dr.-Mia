<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Seurat: R Toolkit for Single Cell Genomics

**Source:** [satijalab/seurat](https://github.com/satijalab/seurat)
**Local Repository:** `./repo`

## Overview
Seurat is the industry-standard R toolkit for single-cell genomics. It provides a comprehensive workflow for QC, analysis, and exploration of scRNA-seq data.

## Installation
To install from the local repository (requires R):
```r
setwd("./repo")
remotes::install_local(".")
```

## Basic Usage (R)
```r
library(Seurat)

# Load Data
pbmc.data <- Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)

# Dimensionality Reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Plotting
DimPlot(pbmc, reduction = "umap")
```

## Tutorials & Vignettes
Seurat includes detailed R Markdown vignettes within its repository.
- **Location:** `./repo/vignettes`
- **Contents:** Rmd files covering standard workflows, integration, multimodal analysis, and spatial genomics.
- **How to use:** Open these `.Rmd` files in RStudio to step through the analysis interactively.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->