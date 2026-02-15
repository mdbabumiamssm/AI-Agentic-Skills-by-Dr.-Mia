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

# Seurat Data: Dataset Manager

**Source:** [satijalab/seurat-data](https://github.com/satijalab/seurat-data)
**Local Repository:** `./repo`

## Overview
SeuratData provides a mechanism to easily install and load datasets for testing and demonstration purposes in Seurat.

## Usage
```r
library(SeuratData)

# Install a dataset
InstallData("pbmc3k")

# Load it
data("pbmc3k")
pbmc3k
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->