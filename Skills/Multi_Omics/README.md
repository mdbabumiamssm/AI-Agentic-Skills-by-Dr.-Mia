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

# Multi-Omics Integration Skills

Skills for integrating and analyzing data across multiple omics modalities.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **multi-omics-integration** | Cross-modality analysis | Factor analysis, network fusion |

## Key Capabilities

- **Factor Analysis** - MOFA, iCluster, SNF
- **Network Integration** - Multi-layer networks
- **Data Fusion** - Joint dimensionality reduction
- **Correlation Analysis** - Cross-omics correlations
- **Pathway Integration** - Multi-omics pathway enrichment

## Key Tools

- **MOFA2** - Multi-Omics Factor Analysis
- **SNF** - Similarity Network Fusion
- **mixOmics** - Multi-omics integration in R
- **omicade4** - Multiple co-inertia analysis
- **PARADIGM** - Pathway-based integration

## Example MOFA Analysis

```python
from mofapy2.run.entry_point import entry_point

# Initialize MOFA
ent = entry_point()

# Set data (multiple views)
ent.set_data_matrix([
    [expression_matrix],  # View 1: RNA-seq
    [methylation_matrix], # View 2: Methylation
    [protein_matrix]      # View 3: Proteomics
], views_names=['RNA', 'Methylation', 'Protein'])

# Run model
ent.build()
ent.run()

# Extract factors
factors = ent.model.getExpectations()['Z']
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->