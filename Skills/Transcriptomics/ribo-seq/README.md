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

# ribo-seq

## Overview

Analyze ribosome profiling (Ribo-seq) data to study translation at single-codon resolution, including periodicity QC, ORF detection, and translation efficiency calculation.

**Tool type:** mixed | **Primary tools:** Plastid, RiboCode, ORFik, riborex

## Skills

| Skill | Description |
|-------|-------------|
| riboseq-preprocessing | Size selection, rRNA removal, and alignment for ribosome footprints |
| ribosome-periodicity | Validate 3-nucleotide periodicity and calculate P-site offsets |
| orf-detection | Detect and quantify translated ORFs with RiboCode and ORFquant |
| translation-efficiency | Calculate translation efficiency from Ribo-seq and RNA-seq |
| ribosome-stalling | Detect ribosome pausing and stalling at specific codons |

## Example Prompts

- "Check ribosome footprint length distribution"
- "Verify 3-nucleotide periodicity at start codons"
- "Find actively translated uORFs"
- "Calculate translation efficiency for my genes"
- "Detect ribosome stalling sites"
- "Generate metagene plots around start codons"
- "Quantify ORF-level translation with ORFquant"
- "Compare ORF expression across conditions"

## Requirements

```bash
# Plastid (Python)
pip install plastid

# RiboCode
pip install RiboCode

# Differential TE and ORF quantification (R)
BiocManager::install(c('riborex', 'Ribo-seQC', 'ORFik'))

# Alignment
conda install -c bioconda star bowtie2 sortmerna
```

## Related Skills

- **read-alignment** - Alignment to transcriptome
- **rna-quantification** - RNA-seq quantification for TE
- **differential-expression** - Compare TE between conditions


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->