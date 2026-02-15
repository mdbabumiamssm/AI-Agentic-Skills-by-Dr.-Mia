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

# Seurat Wrappers: Community Extensions

**Source:** [satijalab/seurat-wrappers](https://github.com/satijalab/seurat-wrappers)
**Local Repository:** `./repo`

## Overview
Seurat Wrappers is a collection of wrapper functions that allow you to run external community tools (like Monocle, Velocyto, LIGER, Harmony) directly on Seurat objects.

## Examples
- **Monocle3:** Trajectory inference.
- **Harmony:** Batch correction.
- **LIGER:** NMF-based integration.

## Usage Example (Harmony)
```r
library(Seurat)
library(SeuratWrappers)
library(harmony)

# Run Harmony integration via wrapper
pbmc <- RunHarmony(pbmc, group.by.vars = "batch")
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->