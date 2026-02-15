<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Azimuth: Reference Mapping

**Source:** [satijalab/azimuth](https://github.com/satijalab/azimuth)
**Local Repository:** `./repo`

## Overview
Azimuth is a Seurat-based tool for mapping query datasets to high-quality reference atlases. It automates annotation by projecting your data onto a pre-computed reference, transferring cell type labels and embeddings.

## Key Features
- **Fast Mapping:** Uses PCA projection for speed.
- **Reference Atlases:** Access to Human PBMC, Cortex, Lung, etc.
- **Web App:** Can be run as a Shiny app or via R command line.

## Basic Usage (R)
```r
library(Seurat)
library(Azimuth)

# Run Azimuth on a query dataset using a built-in reference
pbmc.query <- RunAzimuth(pbmc.query, reference = "pbmcref")
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->