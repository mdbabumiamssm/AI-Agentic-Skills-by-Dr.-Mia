# Imaging Analysis Skills

Skills for biological imaging data analysis including mass cytometry imaging.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **imaging-mass-cytometry** | IMC analysis | Preprocessing, segmentation, spatial stats |

## Key Capabilities

- **Image Preprocessing** - Normalization, background correction
- **Cell Segmentation** - Deep learning-based segmentation
- **Spatial Analysis** - Neighborhood analysis, spatial statistics
- **Phenotyping** - Cell type identification
- **Multiplexed Imaging** - Multi-marker analysis

## Key Tools

- **steinbock** - IMC preprocessing pipeline
- **Cellpose** - Deep learning segmentation
- **squidpy** - Spatial data analysis
- **napari** - Multi-dimensional image viewer
- **scikit-image** - Image processing

## Example IMC Workflow

```python
import scanpy as sc
import squidpy as sq

# Load spatial data
adata = sc.read_h5ad("imc_data.h5ad")

# Neighborhood enrichment analysis
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
sq.pl.nhood_enrichment(adata, cluster_key="cell_type")

# Spatial autocorrelation
sq.gr.spatial_autocorr(adata, mode="moran")

# Co-occurrence analysis
sq.gr.co_occurrence(adata, cluster_key="cell_type")
sq.pl.co_occurrence(adata, cluster_key="cell_type")
```

## Cell Segmentation

```python
from cellpose import models

# Initialize model
model = models.Cellpose(model_type='cyto2')

# Segment cells
masks, flows, styles, diams = model.eval(
    images,
    diameter=30,
    channels=[0, 0],
    flow_threshold=0.4,
    cellprob_threshold=0.0
)
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*
