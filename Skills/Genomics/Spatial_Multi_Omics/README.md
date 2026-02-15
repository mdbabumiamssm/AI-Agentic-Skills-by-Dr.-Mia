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

# Spatial Multi-Omics Agent

**ID:** `biomedical.genomics.spatial_multi_omics`
**Version:** 1.0.0
**Status:** Production
**Category:** Genomics / Spatial Biology / Multi-Omics

---

## Overview

The **Spatial Multi-Omics Agent** enables the analysis of spatially resolved multi-omics data, combining gene expression, protein abundance, chromatin accessibility, and metabolite information with precise tissue location. This agent integrates cutting-edge spatial technologies including MERFISH, Stereo-seq, Visium, CosMx, Xenium, and CODEX to reveal the spatial organization of cellular ecosystems.

Spatial multi-omics is revolutionizing our understanding of tissue architecture, tumor microenvironments, and cell-cell interactions by preserving the critical spatial context lost in dissociated single-cell methods.

---

## Key Capabilities

### 1. Supported Spatial Technologies

| Technology | Company | Resolution | Plex | Modality |
|------------|---------|------------|------|----------|
| **MERFISH** | Vizgen | Subcellular | 500+ genes | RNA |
| **Stereo-seq** | STOmics/BGI | 500nm bins | Whole transcriptome | RNA |
| **Visium HD** | 10x Genomics | 2μm | Whole transcriptome | RNA |
| **Xenium** | 10x Genomics | Subcellular | 5000 genes | RNA |
| **CosMx** | NanoString | Subcellular | 6000+ genes | RNA + Protein |
| **CODEX** | Akoya | Single-cell | 100+ markers | Protein |
| **Stereo-CITE** | STOmics | 500nm | RNA + Protein | Multi-modal |

### 2. Spatial Analysis Capabilities

| Analysis | Description | Tools |
|----------|-------------|-------|
| **Cell Segmentation** | AI-powered cell boundary detection | Cellpose, StarDist, Baysor |
| **Cell Type Mapping** | Spatial cell type annotation | RCTD, Cell2location, Tangram |
| **Spatial Domains** | Tissue region identification | SpaGCN, STAGATE, GraphST |
| **Cell-Cell Interactions** | Ligand-receptor analysis | CellChat, NicheNet, COMMOT |
| **Spatial Trajectories** | Developmental gradients | PAGA, stLearn |
| **Niche Analysis** | Cellular neighborhood composition | Squidpy, CytoSPACE |

### 3. AI-Powered Image Analysis

| Method | Application | Description |
|--------|-------------|-------------|
| **STARVUE** | Vizgen | AI-driven spatial image analysis |
| **Deep learning segmentation** | Cell boundaries | Cellpose 2.0, 3.0 models |
| **Spatial transformers** | Cross-modal integration | SpaFormer, CellPLM spatial |
| **Graph neural networks** | Neighborhood analysis | Spatial GCN, GAT |

---

## Usage

### Example Prompt

```text
Analyze this Xenium breast cancer spatial transcriptomics dataset to characterize the tumor microenvironment.

Data:
- Xenium data: breast_tumor_xenium/
- 5000-plex gene panel
- FFPE tissue section

Tasks:
1. Perform cell segmentation and quality control
2. Annotate cell types spatially
3. Identify spatial domains/niches
4. Analyze cell-cell interactions
5. Characterize the tumor-immune interface
6. Identify prognostic spatial patterns
```

### Expected Output

```
## Spatial Transcriptomics Analysis Report

### Dataset Summary
- **Technology:** 10x Genomics Xenium
- **Sample:** Breast cancer FFPE
- **Tissue area:** 12.4 mm²
- **Genes profiled:** 5,000
- **Transcripts detected:** 28.7 million

### Cell Segmentation Results

| Metric | Value |
|--------|-------|
| Total cells | 234,567 |
| Median transcripts/cell | 122 |
| Median genes/cell | 87 |
| Segmentation method | Cellpose + Xenium nuclei |

### Spatial Cell Type Annotation

| Cell Type | Count | % | Spatial Distribution |
|-----------|-------|---|---------------------|
| Tumor cells | 89,234 | 38.0% | Core and invasive front |
| CD8+ T cells | 23,456 | 10.0% | Tumor margin, tertiary lymphoid |
| CD4+ T cells | 18,234 | 7.8% | T cell zones |
| B cells | 12,345 | 5.3% | Tertiary lymphoid structures |
| Macrophages (M1) | 15,678 | 6.7% | Tumor-stroma interface |
| Macrophages (M2) | 21,234 | 9.1% | Within tumor core |
| Fibroblasts (CAF) | 28,456 | 12.1% | Stromal regions |
| Endothelial | 8,934 | 3.8% | Vascular structures |
| Other | 16,996 | 7.2% | - |

### Spatial Domain Analysis

**6 Distinct Spatial Domains Identified:**

| Domain | Characteristics | Cells | Area |
|--------|-----------------|-------|------|
| **Tumor Core** | High proliferation, hypoxic | 67,234 | 3.2 mm² |
| **Invasive Front** | EMT signature, matrix remodeling | 34,567 | 1.8 mm² |
| **Immune Hot** | T cell rich, TLS adjacent | 28,456 | 1.4 mm² |
| **Immune Cold** | T cell excluded, immunosuppressive | 45,678 | 2.1 mm² |
| **Stroma** | CAF-rich, ECM deposition | 38,234 | 2.6 mm² |
| **Normal Adjacent** | Normal epithelium, low immune | 20,398 | 1.3 mm² |

### Tertiary Lymphoid Structures (TLS)

| Metric | Value |
|--------|-------|
| TLS detected | 12 structures |
| Mean TLS area | 0.08 mm² |
| TLS maturity | 8 mature, 4 immature |
| B cell follicles | 8 present |
| T cell zones | All 12 |

**TLS Spatial Association:**
- TLS within 200μm of tumor: 9/12 (75%)
- Associated with better prognosis (HR=0.45, p<0.001)

### Cell-Cell Interaction Analysis

#### Top Ligand-Receptor Pairs (Tumor-Immune Interface)

| Rank | Ligand | Receptor | Sender | Receiver | Score |
|------|--------|----------|--------|----------|-------|
| 1 | PD-L1 | PD-1 | Tumor | CD8+ T | 0.89 |
| 2 | CXCL9 | CXCR3 | Macrophage | CD8+ T | 0.84 |
| 3 | CCL19 | CCR7 | FDC | B cell | 0.81 |
| 4 | SPP1 | CD44 | Macrophage | Tumor | 0.78 |
| 5 | COL1A1 | ITGA2 | CAF | Tumor | 0.76 |

### Tumor-Immune Interface Characterization

#### Spatial Immune Phenotype Classification

| Phenotype | Criteria | % of Tumor |
|-----------|----------|------------|
| **Inflamed** | CD8+ within tumor | 32% |
| **Excluded** | CD8+ at margin only | 48% |
| **Desert** | Minimal CD8+ anywhere | 20% |

#### Immune Exclusion Zone Analysis
- Mean exclusion zone width: 145 μm
- CAF density in zone: 3.2x higher than stroma
- SPP1+ macrophages enriched at interface

### Prognostic Spatial Patterns

| Pattern | Hazard Ratio | 95% CI | P-value |
|---------|--------------|--------|---------|
| TLS density (high vs low) | 0.45 | 0.28-0.72 | 0.001 |
| CD8+ tumor infiltration | 0.52 | 0.34-0.79 | 0.002 |
| CAF barrier thickness | 2.34 | 1.56-3.51 | <0.001 |
| M2/M1 ratio (high) | 1.89 | 1.23-2.90 | 0.003 |

### Spatial Biomarker Summary

**Favorable Spatial Features:**
- TLS presence and maturity
- CD8+ T cell infiltration into tumor
- M1-polarized macrophages at interface
- Thin/absent CAF barrier

**Unfavorable Spatial Features:**
- Immune exclusion phenotype
- Thick CAF barrier (>200 μm)
- M2 macrophage dominance in core
- Absence of TLS

### Visualization Outputs

Generated files:
1. `cell_type_spatial_map.png` - Cell type overlay
2. `spatial_domains.png` - Domain segmentation
3. `ligand_receptor_network.html` - Interactive L-R
4. `tls_analysis.png` - TLS characterization
5. `immune_phenotype_map.png` - Inflamed/excluded/desert
```

### LLM Agent Integration

```python
@tool
def segment_spatial_cells(
    spatial_data_dir: str,
    method: str = "cellpose",
    model: str = "cyto3",
    use_nuclei: bool = True
) -> str:
    """
    Segments cells from spatial transcriptomics images.

    Args:
        spatial_data_dir: Path to spatial data directory
        method: Segmentation method (cellpose, stardist, baysor)
        model: Model for segmentation
        use_nuclei: Use nuclei channel for guidance

    Returns:
        Segmentation results with cell boundaries
    """
    pass


@tool
def annotate_spatial_cell_types(
    spatial_data: str,
    reference: str = "human_cell_atlas",
    method: str = "cell2location",
    return_confidence: bool = True
) -> str:
    """
    Annotates cell types in spatial context.

    Args:
        spatial_data: Path to processed spatial data
        reference: Reference scRNA-seq atlas
        method: Annotation method (cell2location, rctd, tangram)
        return_confidence: Include confidence scores

    Returns:
        Spatial cell type annotations
    """
    pass


@tool
def identify_spatial_domains(
    spatial_data: str,
    n_domains: int = None,
    method: str = "stagate",
    resolution: float = 1.0
) -> str:
    """
    Identifies spatial domains/tissue regions.

    Args:
        spatial_data: Path to spatial data
        n_domains: Number of domains (None = auto)
        method: Clustering method (stagate, spagcn, graphst)
        resolution: Clustering resolution

    Returns:
        Domain assignments with characterization
    """
    pass


@tool
def analyze_cell_interactions(
    spatial_data: str,
    method: str = "cellchat",
    max_distance: float = 200.0,
    ligand_receptor_db: str = "cellchatdb"
) -> str:
    """
    Analyzes spatially-constrained cell-cell interactions.

    Args:
        spatial_data: Path to spatial data
        method: Interaction analysis (cellchat, nichenet, commot)
        max_distance: Maximum interaction distance (μm)
        ligand_receptor_db: L-R database to use

    Returns:
        Significant cell-cell interactions with spatial context
    """
    pass


@tool
def characterize_tumor_microenvironment(
    spatial_data: str,
    tumor_annotation: str,
    analyses: list[str] = ["immune_phenotype", "tls", "interface"]
) -> str:
    """
    Comprehensive TME characterization from spatial data.

    Args:
        spatial_data: Path to spatial data
        tumor_annotation: Column with tumor cell labels
        analyses: TME analyses to perform

    Returns:
        TME report with spatial patterns
    """
    pass
```

---

## Prerequisites

### Required Packages

| Package | Version | Purpose |
|---------|---------|---------|
| **squidpy** | >=1.3 | Spatial analysis |
| **scanpy** | >=1.9 | Single-cell analysis |
| **cellpose** | >=3.0 | Cell segmentation |
| **cell2location** | >=0.1 | Cell type deconvolution |
| **stlearn** | >=0.4 | Spatial trajectories |

### Dependencies

```
squidpy>=1.3.0
scanpy>=1.9.0
cellpose>=3.0.0
cell2location>=0.1.0
spatialdata>=0.1.0
napari>=0.4
torch>=2.0
pandas>=2.0
numpy>=1.24
```

---

## Methodology

### Spatial Analysis Pipeline

```
Raw Spatial Data (images + transcripts)
    ↓
Image Processing
├── Tissue detection
├── Nuclear staining
└── Transcript localization
    ↓
Cell Segmentation
├── Cellpose/StarDist
├── Nuclear expansion
└── Transcript assignment
    ↓
Quality Control
├── Cell size filtering
├── Transcript QC
└── Doublet detection
    ↓
Cell Type Annotation
├── Reference mapping
├── Marker-based
└── Deconvolution
    ↓
Spatial Analysis
├── Domain identification
├── Neighborhood analysis
├── Cell-cell interactions
└── Spatial statistics
    ↓
Visualization & Reporting
```

---

## Clinical Applications

### Cancer Research
- Tumor microenvironment characterization
- Immunotherapy response prediction
- Spatial biomarker discovery
- Metastasis mechanism studies

### Tissue Biology
- Organ atlas construction
- Development studies
- Disease tissue mapping
- Cell fate in situ

### Drug Development
- Target validation in tissue context
- Drug distribution studies
- Mechanism of action in situ

---

## Related Skills

- **Single-Cell Analysis Agent:** Dissociated scRNA-seq
- **Tumor Microenvironment Agent:** TME profiling
- **Digital Pathology Agent:** H&E analysis
- **Multi-Omics Integration:** Cross-modal analysis

---

## References

- **Chen et al. (2022):** "Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays." *Cell*
- **He et al. (2022):** "High-plex imaging of RNA and proteins at subcellular resolution in fixed tissue by spatial molecular imaging." *Nature Biotechnology*
- [Squidpy Documentation](https://squidpy.readthedocs.io/)
- [Vizgen MERFISH](https://vizgen.com/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->