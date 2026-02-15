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

# hi-c-analysis

## Overview

Analyze Hi-C and chromosome conformation capture data using cooler, cooltools, pairtools, and HiCExplorer.

**Tool type:** mixed | **Primary tools:** cooler, cooltools, pairtools, HiCExplorer

## Skills

| Skill | Description |
|-------|-------------|
| hic-data-io | Load, convert, and manipulate Hi-C matrices in cooler format |
| contact-pairs | Process Hi-C read pairs with pairtools |
| matrix-operations | Balance, normalize, and transform contact matrices |
| compartment-analysis | Detect A/B compartments from Hi-C data |
| tad-detection | Call topologically associating domains (TADs) |
| loop-calling | Detect chromatin loops and point interactions |
| hic-visualization | Visualize contact matrices and genomic features |
| hic-differential | Compare Hi-C between conditions, differential contacts |

## Example Prompts

- "Load my Hi-C contact matrix"
- "Convert .hic to cooler format"
- "Process Hi-C read pairs with pairtools"
- "Balance my Hi-C matrix"
- "Call compartments from my Hi-C data"
- "Detect TADs in my contact matrix"
- "Find chromatin loops"
- "Plot a Hi-C contact matrix"
- "Show interactions around this gene"
- "Compare Hi-C matrices between conditions"
- "Find differential contacts between treatment and control"
- "Identify compartment switches between conditions"

## Requirements

```bash
# Python tools
pip install cooler cooltools bioframe matplotlib

# CLI tools (via conda)
conda install -c bioconda pairtools hicexplorer fanc
```

## Related Skills

- **genome-intervals** - Work with genomic coordinates
- **chip-seq** - Integrate with ChIP-seq peaks
- **atac-seq** - Integrate with accessibility data


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->