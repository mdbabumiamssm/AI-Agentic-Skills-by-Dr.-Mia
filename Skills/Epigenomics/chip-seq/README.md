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

# chip-seq

## Overview

ChIP-seq analysis using MACS3 for peak calling, ChIPseeker for annotation, and DiffBind for differential binding. Used to study protein-DNA interactions including transcription factor binding and histone modifications.

**Tool type:** mixed | **Primary tools:** MACS3 (CLI), ChIPseeker (R), DiffBind (R)

## Skills

| Skill | Description |
|-------|-------------|
| peak-calling | Call peaks with MACS3 from ChIP-seq BAM files |
| peak-annotation | Annotate peaks to genes with ChIPseeker |
| differential-binding | Differential binding analysis with DiffBind |
| chipseq-visualization | Visualize ChIP-seq data and peaks |
| motif-analysis | De novo and known motif enrichment with HOMER/MEME |
| chipseq-qc | Quality metrics: FRiP, NSC/RSC, IDR for replicate concordance |
| super-enhancers | Super-enhancer identification with ROSE |

## Example Prompts

- "Call peaks from my ChIP-seq BAM with MACS3"
- "Run MACS3 with input control"
- "Call broad peaks for H3K27me3"
- "Annotate my peaks with nearby genes"
- "Find peaks in promoters vs enhancers"
- "What genes have peaks in their promoter?"
- "Find differential peaks between treatment and control"
- "Compare binding between two conditions with DiffBind"
- "Identify lost and gained peaks"
- "Create a heatmap of peak signal"
- "Plot peak distribution around TSS"
- "Visualize peaks in a genomic region"
- "Find enriched motifs in my peaks with HOMER"
- "Run de novo motif discovery"
- "What transcription factors bind my peaks?"
- "Calculate FRiP for my ChIP-seq experiment"
- "Run IDR on my replicates"
- "Check cross-correlation (NSC/RSC) for my ChIP"

## Requirements

```bash
# MACS3
conda install -c bioconda macs3

# HOMER (motif analysis)
conda install -c bioconda homer

# MEME Suite (alternative motif analysis)
conda install -c bioconda meme

# QC tools
conda install -c bioconda phantompeakqualtools idr deeptools
```

```r
BiocManager::install(c('ChIPseeker', 'DiffBind', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

## Related Skills

- **alignment-files** - BAM preparation before peak calling
- **pathway-analysis** - Functional enrichment of peak-associated genes
- **genome-intervals** - BED file operations for peak regions
- **methylation-analysis** - Combine with DNA methylation data


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->