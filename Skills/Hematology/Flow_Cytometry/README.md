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

# flow-cytometry

## Overview

Flow and mass cytometry analysis from FCS files to differential cell populations.

**Tool type:** r | **Primary tools:** flowCore, CATALYST, CytoML

## Skills

| Skill | Description |
|-------|-------------|
| fcs-handling | Read and write FCS files, access channel data |
| compensation-transformation | Spillover compensation and data transformation |
| gating-analysis | Manual and automated gating strategies |
| clustering-phenotyping | Unsupervised clustering with FlowSOM and Phenograph |
| differential-analysis | Differential abundance and state analysis |
| doublet-detection | Remove doublet events using FSC-A vs FSC-H |
| bead-normalization | CyTOF EQ bead normalization and drift correction |
| cytometry-qc | Comprehensive QC for flow and mass cytometry |

## Example Prompts

- "Load my FCS files and apply compensation"
- "Create a gating strategy for T cell subsets"
- "Cluster my CyTOF data with FlowSOM"
- "Find differentially abundant populations between groups"
- "Perform dimension reduction with UMAP"

## Requirements

```r
# R/Bioconductor
BiocManager::install(c("flowCore", "flowWorkspace", "CytoML", "CATALYST", "FlowSOM", "ggcyto"))
```

## Related Skills

- **single-cell** - Similar analysis concepts
- **imaging-mass-cytometry** - Related spatial proteomics
- **multi-omics-integration** - Integrate with other data types


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->