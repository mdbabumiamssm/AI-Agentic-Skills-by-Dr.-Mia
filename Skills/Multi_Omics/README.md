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
