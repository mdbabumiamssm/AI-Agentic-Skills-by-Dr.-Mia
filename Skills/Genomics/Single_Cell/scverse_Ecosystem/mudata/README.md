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

# MuData: Multimodal Data Containers

**Source:** [scverse/mudata](https://github.com/scverse/mudata)
**Local Repository:** `./repo`

## Overview
MuData provides a container for multimodal omics data (e.g., CITE-seq, Multiome ATAC+RNA). It acts as a wrapper around multiple `AnnData` objects, keeping them synchronized.

## Structure
A `MuData` object contains a dictionary of `AnnData` objects (modalities).
- `mdata.mod['rna']`: AnnData for RNA.
- `mdata.mod['atac']`: AnnData for ATAC.

## Basic Usage
```python
import mudata as md
from anndata import AnnData
import numpy as np

# Create multimodal data
rna = AnnData(np.random.normal(size=(100, 1000)))
atac = AnnData(np.random.normal(size=(100, 500)))

mdata = md.MuData({"rna": rna, "atac": atac})

# Update logic
mdata.update() # Syncs obs_names across modalities
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->