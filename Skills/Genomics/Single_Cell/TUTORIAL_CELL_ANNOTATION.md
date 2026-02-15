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

# Tutorial: Single-Cell Cell Type Annotation

This tutorial guides you through using the `UniversalAnnotator` skill to identify cell types in single-cell RNA-sequencing (scRNA-seq) data.

## Prerequisites

Ensure you have the necessary libraries installed. While the script has a mock fallback, real analysis requires:

```bash
pip install scanpy pandas numpy celltypist
```

## Scenario

You have a gene expression matrix (cells x genes). You want to label the cells based on their gene expression profiles. We will use a hybrid approach:
1.  **Marker Scoring**: Fast, rule-based.
2.  **CellTypist**: Deep Learning based, high accuracy for immune cells.

## Step 1: Basic Usage (Marker Based)

Create a Python script `annotate_markers.py`:

```python
import scanpy as sc
import numpy as np
from Skills.Genomics.Single_Cell.Cell_Type_Annotation.RNA.universal_annotator import UniversalAnnotator

# 1. Load your data (or generate mock data)
adata = sc.AnnData(np.random.rand(100, 2000))
# Assign some gene names matching our markers
adata.var_names = [f"Gene_{i}" for i in range(2000)]
# Force some genes to exist for the demo
idx_map = {0: 'CD3D', 1: 'CD14', 2: 'MS4A1'}
new_names = list(adata.var_names)
for i, name in idx_map.items():
    new_names[i] = name
adata.var_names = new_names

# 2. Initialize Annotator
annotator = UniversalAnnotator(adata)

# 3. Define Markers
markers = {
    'T-Cells': ['CD3D'],
    'Monocytes': ['CD14'],
    'B-Cells': ['MS4A1']
}

# 4. Run Annotation
annotator.annotate_marker_based(markers)

# 5. View Results
print(adata.obs[['predicted_cell_type']].value_counts())
```

Run it:
```bash
python3 annotate_markers.py
```

## Step 2: Advanced Usage (CellTypist)

If you have `celltypist` installed, you can leverage pre-trained models.

```python
# ... load adata ...

# Run CellTypist
annotator.annotate_with_celltypist(model_name='Immune_All_Low.pkl')

# Check overlap between Marker predictions and CellTypist predictions
pd.crosstab(adata.obs['predicted_cell_type'], adata.obs['celltypist_prediction'])
```

## Step 3: LLM-Assisted Discovery

For clusters that don't match known markers or CellTypist models, use the LLM workflow.

```python
# 1. Cluster the data first
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added='leiden')

# 2. Extract markers and print prompts
annotator.annotate_with_llm(cluster_col='leiden', marker_num=5)

# 3. Copy-paste the output prompts into ChatGPT/Claude or automate via API.
```

## Summary

*   **Marker-based**: Good for verifying known biology.
*   **CellTypist**: Best for standardized immune phenotyping.
*   **LLM**: Best for exploring unknown clusters.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->