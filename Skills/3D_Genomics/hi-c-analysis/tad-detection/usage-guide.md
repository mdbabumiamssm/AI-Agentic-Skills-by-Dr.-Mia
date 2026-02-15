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

# TAD Detection - Usage Guide

## Overview

This skill covers calling topologically associating domains (TADs) from Hi-C data using insulation score methods with cooltools and HiCExplorer.

## Prerequisites

```bash
pip install cooler cooltools bioframe matplotlib
# For HiCExplorer:
conda install -c bioconda hicexplorer
```

## Quick Start

Tell your AI agent what you want to do:

- "Call TADs from my Hi-C data"
- "Find TAD boundaries using insulation score"
- "Compute domain boundaries at 10kb resolution"
- "Compare TAD structures between conditions"

## Example Prompts

### TAD Calling
> "Detect TADs from this cooler file"

> "Compute insulation score and find boundaries"

### Analysis
> "What is the average TAD size?"

> "Compare TAD boundaries between samples"

### Export
> "Save TAD boundaries as BED file"

## What the Agent Will Do

1. Load cooler at appropriate resolution (10-25kb)
2. Compute insulation score
3. Identify local minima as boundaries
4. Define TADs between boundaries
5. Return TAD intervals and statistics

## Tips

- **Resolution** - Use 10-25kb for TAD calling
- **Window size** - 100-500kb typical; larger = fewer boundaries
- **Boundary strength** - Filter weak boundaries
- **Multiple windows** - Reveals hierarchical structure


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->