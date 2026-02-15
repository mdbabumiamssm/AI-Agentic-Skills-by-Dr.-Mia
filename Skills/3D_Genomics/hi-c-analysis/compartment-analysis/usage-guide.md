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

# Compartment Analysis - Usage Guide

## Overview

This skill covers detecting A/B compartments from Hi-C data using eigenvector decomposition with cooltools.

## Prerequisites

```bash
pip install cooler cooltools bioframe matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Call compartments from my Hi-C data"
- "Compute A/B compartments"

## Example Prompts

### Compartment Calling
> "Detect compartments from this cooler file"

> "Compute the first eigenvector for compartment analysis"

### Visualization
> "Plot a saddle plot for compartment strength"

> "Show the compartment track for chr1"

### Comparison
> "Compare compartments between treatment and control"

## What the Agent Will Do

1. Load cooler at appropriate resolution (50-100kb)
2. Compute expected values
3. Compute eigenvector decomposition
4. Assign A/B compartments based on E1 sign
5. Optionally compute saddle plot

## Tips

- **Resolution** - Use 50-100kb for compartment analysis
- **GC phasing** - Use GC content to correctly orient A/B
- **E1 sign** - Positive typically = A (active), negative = B
- **Saddle plot** - Shows compartmentalization strength


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->