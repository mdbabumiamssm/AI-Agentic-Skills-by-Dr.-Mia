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

# Epigenomics Skills

Skills for epigenetic data analysis including ChIP-seq, ATAC-seq, and methylation.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **chip-seq** | ChIP-seq analysis | Peak calling, annotation, differential binding |
| **atac-seq** | ATAC-seq analysis | Chromatin accessibility, footprinting |
| **methylation-analysis** | DNA methylation | Bisulfite alignment, DMR detection |
| **epitranscriptomics** | RNA modifications | m6A analysis, MeRIP-seq |
| **clip-seq** | CLIP-seq analysis | Protein-RNA interactions |

## Key Tools

- **MACS3** - Peak calling
- **deepTools** - Signal visualization
- **ArchR** - Single-cell ATAC-seq
- **Bismark** - Bisulfite alignment
- **DSS/DMRcate** - DMR detection

## Example Workflow

```bash
# ChIP-seq peak calling
macs3 callpeak -t treatment.bam -c control.bam \
    -f BAM -g hs -n experiment --outdir peaks/

# ATAC-seq analysis
macs3 callpeak -t atac.bam -f BAM -g hs -n atac \
    --nomodel --shift -100 --extsize 200
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->