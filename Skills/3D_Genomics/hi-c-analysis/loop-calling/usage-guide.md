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

# Loop Calling - Usage Guide

## Overview

This skill covers detecting chromatin loops and point interactions from Hi-C data using cooltools, chromosight, and HiCExplorer.

## Prerequisites

```bash
pip install cooler cooltools bioframe matplotlib
# For chromosight:
pip install chromosight
# For HiCExplorer:
conda install -c bioconda hicexplorer
```

## Quick Start

Tell your AI agent what you want to do:
- "Call loops from my Hi-C data"
- "Find chromatin loops"

## Example Prompts

### Loop Calling
> "Detect loops from this cooler file"

> "Find point interactions in my Hi-C data"

### Analysis
> "What is the average loop size?"

> "How many loops have CTCF at both anchors?"

### APA
> "Compute aggregate peak analysis"

## What the Agent Will Do

1. Load cooler at 5-10kb resolution
2. Compute expected values
3. Call loops using dot detection
4. Filter by score/significance
5. Return loop coordinates

## Tips

- **Resolution** - Use 5-10kb for loop calling
- **Max distance** - Most loops are <2Mb
- **CTCF** - Validate loops against CTCF peaks
- **APA** - Stack loops to confirm enrichment


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->