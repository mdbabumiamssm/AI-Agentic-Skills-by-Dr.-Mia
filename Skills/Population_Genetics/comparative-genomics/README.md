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

# comparative-genomics

## Overview

Comparative genomics including synteny analysis, positive selection detection, ancestral sequence reconstruction, ortholog inference, and horizontal gene transfer detection.

**Tool type:** mixed | **Primary tools:** MCScanX, PAML, OrthoFinder, HyPhy

## Skills

| Skill | Description |
|-------|-------------|
| synteny-analysis | Genome collinearity with MCScanX, SyRI, JCVI |
| positive-selection | dN/dS tests with PAML codeml, HyPhy |
| ancestral-reconstruction | Ancestral sequences with PAML, IQ-TREE |
| ortholog-inference | Orthogroups with OrthoFinder, ProteinOrtho |
| hgt-detection | Horizontal gene transfer with HGTector |

## Example Prompts

- "Find syntenic blocks between human and mouse genomes"
- "Test for positive selection on this gene across mammals"
- "Reconstruct the ancestral sequence at this node"
- "Find orthologs of BRCA1 across vertebrates"
- "Detect horizontally transferred genes in this bacterial genome"
- "Calculate dN/dS ratios for this gene family"

## Requirements

```bash
# OrthoFinder
conda install -c bioconda orthofinder

# PAML
conda install -c bioconda paml

# HyPhy (modern alternative to PAML)
conda install -c bioconda hyphy

# MCScanX
git clone https://github.com/wyp1125/MCScanX
cd MCScanX && make
```

## Related Skills

- **phylogenetics** - Tree building for selection analysis
- **alignment** - MSA preparation for codon analysis
- **population-genetics** - Within-species selection statistics


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->