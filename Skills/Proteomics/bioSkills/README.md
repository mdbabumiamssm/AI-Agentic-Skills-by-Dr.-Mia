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

# proteomics

## Overview

Mass spectrometry-based proteomics analysis from raw data to differential abundance.

**Tool type:** mixed | **Primary tools:** pyOpenMS, MSstats, limma, QFeatures

## Skills

| Skill | Description |
|-------|-------------|
| data-import | Load and parse mass spec data formats (mzML, mzXML, MaxQuant output) |
| peptide-identification | Peptide-spectrum matching and protein identification |
| quantification | Label-free and labeled (TMT/iTRAQ/SILAC) protein quantification |
| proteomics-qc | Quality control, missing value analysis, batch effect detection |
| differential-abundance | Statistical testing for differentially abundant proteins |
| ptm-analysis | Post-translational modification identification and localization |
| protein-inference | Protein grouping, inference, and FDR control |
| dia-analysis | Data-independent acquisition analysis with DIA-NN |
| spectral-libraries | Build and use spectral libraries (empirical and predicted) |

## Example Prompts

- "Load MaxQuant proteinGroups.txt and filter contaminants"
- "Perform label-free quantification from mzML files"
- "Find differentially abundant proteins between treatment and control"
- "Identify phosphorylation sites from enriched samples"
- "Calculate protein-level FDR using target-decoy approach"

## Requirements

```bash
# Python
pip install pyopenms pandas numpy scipy

# R
install.packages(c("MSstats", "limma"))
BiocManager::install(c("QFeatures", "proDA", "MSnbase"))
```

## Related Skills

- **differential-expression** - Similar statistical approaches for RNA-seq
- **pathway-analysis** - Functional enrichment of protein lists
- **data-visualization** - Volcano plots, heatmaps for proteomics results


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->