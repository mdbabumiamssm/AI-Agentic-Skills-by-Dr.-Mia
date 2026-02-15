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

# Single-Cell Clustering - Usage Guide

## Overview

This skill covers dimensionality reduction and clustering for single-cell RNA-seq data using both Seurat (R) and Scanpy (Python). These steps identify cell populations in your data.

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy matplotlib leidenalg
```

**R (Seurat):**
```r
install.packages('Seurat')
# Optional for resolution comparison:
install.packages('clustree')
```

## Quick Start

Ask your AI agent:

> "Run PCA and cluster my single-cell data"

> "Generate a UMAP and color by cluster"

> "Try different clustering resolutions"

## Example Prompts

### Dimensionality Reduction
> "Run PCA and show the elbow plot"

> "How many PCs should I use?"

> "Generate a UMAP embedding"

### Clustering
> "Cluster at resolution 0.5"

> "How many cells are in each cluster?"

> "Try resolutions 0.2, 0.5, and 1.0 and compare"

### Visualization
> "Show clusters on UMAP"

> "Color UMAP by sample and by cluster"

> "Plot expression of CD3D on the UMAP"

## What the Agent Will Do

1. Run PCA on variable genes
2. Determine optimal number of PCs
3. Build neighbor graph
4. Run clustering algorithm (Leiden/Louvain)
5. Generate UMAP embedding
6. Visualize results

## Tips

- **Leiden is preferred** over Louvain (better performance, more consistent)
- **Resolution is key** - try multiple values and use biological knowledge
- **10-50 PCs** typically sufficient - check elbow plot
- **n_neighbors affects smoothness** - higher = smoother clusters
- **Save intermediate results** - clustering can take time on large datasets
- **UMAP is for visualization** - don't cluster on UMAP coordinates


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->