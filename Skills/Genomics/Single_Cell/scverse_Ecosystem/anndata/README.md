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

# AnnData: Annotated Data Matrix

**Source:** [scverse/anndata](https://github.com/scverse/anndata)
**Local Repository:** `./repo`

## Overview
AnnData is the core data structure for the scverse ecosystem. It stores the data matrix (`.X`) alongside annotations for observations (cells, `.obs`) and variables (genes, `.var`).

## Structure
- `adata.X`: The main data matrix (numpy array or sparse matrix).
- `adata.obs`: Dataframe of cell annotations (e.g., cell type, batch).
- `adata.var`: Dataframe of gene annotations (e.g., gene name).
- `adata.obsm`: Multi-dimensional observation annotations (e.g., PCA, UMAP coordinates).
- `adata.uns`: Unstructured data (metadata).

## Basic Usage
```python
import anndata as ad
import pandas as pd
import numpy as np

# Create an AnnData object
counts = np.random.randint(0, 100, (100, 1000))
obs = pd.DataFrame(index=[f"cell_{i}" for i in range(100)])
var = pd.DataFrame(index=[f"gene_{i}" for i in range(1000)])

adata = ad.AnnData(X=counts, obs=obs, var=var)

# Access data
print(adata.X)
print(adata.obs_names)
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->