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

# liquid-biopsy

## Overview
Cell-free DNA and circulating tumor DNA analysis for non-invasive cancer detection, tumor fraction estimation, mutation detection, and treatment monitoring from plasma samples.

**Tool type:** mixed | **Primary tools:** ichorCNA, FinaleToolkit, fgbio, VarDict

## Skills
| Skill | Description |
|-------|-------------|
| cfdna-preprocessing | Preprocess cfDNA reads with UMI-aware deduplication |
| fragment-analysis | Analyze fragmentomics patterns for cancer detection |
| tumor-fraction-estimation | Estimate ctDNA fraction from shallow WGS |
| ctdna-mutation-detection | Detect somatic mutations at low VAF |
| methylation-based-detection | Analyze cfDNA methylation for cancer detection |
| longitudinal-monitoring | Track ctDNA dynamics over treatment |

## Example Prompts
- "Preprocess my plasma cfDNA FASTQ files with UMI handling"
- "Analyze fragment size distribution for tumor signal"
- "Estimate tumor fraction from my shallow WGS data"
- "Detect mutations at 0.5% VAF from my targeted panel"
- "Track ctDNA levels across my serial samples"

## Requirements
```bash
# Python
pip install finaletoolkit pysam pandas matplotlib

# R
BiocManager::install('ichorCNA')

# CLI
conda install -c bioconda fgbio vardict
```

## Related Skills

- **variant-calling** - Somatic variant calling principles
- **copy-number** - CNV detection concepts
- **methylation-analysis** - Methylation processing


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->