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

# microbiome

## Overview

16S/ITS amplicon sequencing analysis from raw reads to differential abundance testing.

**Tool type:** mixed | **Primary tools:** DADA2, phyloseq, ALDEx2, QIIME2

## Skills

| Skill | Description |
|-------|-------------|
| amplicon-processing | ASV inference from 16S/ITS amplicon reads with DADA2 |
| taxonomy-assignment | Taxonomic classification with SILVA, GTDB, or UNITE databases |
| diversity-analysis | Alpha and beta diversity metrics with phyloseq |
| differential-abundance | Compositional-aware testing with ALDEx2, ANCOM-BC |
| functional-prediction | Predict metagenome function from 16S with PICRUSt2 |
| qiime2-workflow | QIIME2 CLI-based amplicon analysis workflow |

## Example Prompts

- "Process paired-end 16S reads and infer ASVs with DADA2"
- "Assign taxonomy using SILVA 138 database"
- "Calculate alpha diversity (Shannon, Chao1) and compare groups"
- "Find differentially abundant taxa between treatment and control"
- "Predict KEGG pathways from 16S data with PICRUSt2"

## Requirements

```bash
# R
install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq", "ALDEx2", "ANCOMBC"))

# QIIME2 (conda)
conda install -c qiime2 qiime2

# PICRUSt2 (conda)
conda install -c bioconda picrust2
```

## Related Skills

- **metagenomics** - Shotgun metagenomics (taxonomic and functional)
- **pathway-analysis** - Enrichment analysis of predicted functions
- **data-visualization** - Ordination plots, taxonomic barplots


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->