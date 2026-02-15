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

# epitranscriptomics

## Overview

Analysis of RNA modifications (m6A, m5C) from MeRIP-seq and direct RNA sequencing.

**Tool type:** mixed | **Primary tools:** exomePeak2, MACS3, m6Anet, Guitar

## Skills

| Skill | Description |
|-------|-------------|
| merip-preprocessing | Align MeRIP-seq IP and input samples with QC |
| m6a-peak-calling | Call m6A peaks from MeRIP-seq data |
| m6a-differential | Identify differential m6A methylation between conditions |
| m6anet-analysis | Detect m6A from ONT direct RNA sequencing |
| modification-visualization | Create metagene plots and browser tracks |

## Example Prompts

- "Align my MeRIP-seq IP and input samples"
- "Call m6A peaks from my MeRIP-seq data"
- "Find differential m6A sites between treatment and control"
- "Detect m6A modifications from my Nanopore direct RNA data"
- "Create a metagene plot of m6A distribution around stop codons"

## Requirements

```bash
# R packages
BiocManager::install(c('exomePeak2', 'Guitar', 'GenomicFeatures'))

# Python packages
pip install m6anet ont-fast5-api

# CLI tools
conda install -c bioconda star macs3 samtools
```

## Related Skills

- **chip-seq** - Similar peak-calling concepts
- **rna-quantification** - RNA-seq alignment
- **long-read-sequencing** - Nanopore data processing


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->