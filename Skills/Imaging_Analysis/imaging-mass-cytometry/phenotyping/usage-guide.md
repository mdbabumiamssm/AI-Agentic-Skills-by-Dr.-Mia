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

# Phenotyping - Usage Guide

## Overview
Phenotyping assigns cell type identities based on marker expression, enabling biological interpretation of IMC data.

## Prerequisites
```bash
pip install scanpy anndata scikit-learn
# For FlowSOM-style clustering
pip install minisom
```

## Quick Start
Tell your AI agent what you want to do:
- "Cluster cells and identify cell types from my IMC data"
- "Assign cell phenotypes using marker expression"
- "Create a manual gating strategy for my panel"

## Example Prompts

### Clustering-Based Phenotyping
> "Cluster my IMC cells using Leiden clustering and annotate cell types"

> "Run FlowSOM clustering on my single-cell data and identify populations"

> "Perform PCA and UMAP on my cell data then cluster with Leiden"

### Manual Gating
> "Define CD45+CD3+CD8+ cells as cytotoxic T cells in my data"

> "Create a gating hierarchy for immune cell populations"

### Marker-Based Annotation
> "Identify macrophages as CD45+CD68+ cells"

> "Annotate epithelial cells using E-cadherin and pan-cytokeratin"

### Validation
> "Show me marker expression heatmaps for each cluster"

> "Validate my cell type assignments against known marker profiles"

## What the Agent Will Do
1. Load single-cell expression data from segmented images
2. Transform data (arcsinh with cofactor 5)
3. Select lineage markers for phenotyping
4. Run dimensionality reduction (PCA, UMAP)
5. Perform clustering (Leiden, FlowSOM)
6. Annotate clusters based on marker expression
7. Generate validation plots (heatmaps, dotplots)

## Tips
- Always arcsinh transform data before clustering
- Select only lineage markers for phenotyping (exclude functional markers)
- Common lineage markers: CD45 (immune), CD3/CD4/CD8 (T cells), CD20 (B cells), CD68 (macrophages)
- FlowSOM produces more reproducible clusters than Leiden for cytometry data
- Validate assignments by checking marker expression in each cluster
- Report confidence of ambiguous cell type assignments


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->