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

# Hi-C Differential Analysis - Usage Guide

## Overview

This skill covers comparing Hi-C contact matrices between conditions to identify differential chromatin interactions, compartment switches, and boundary changes.

## Prerequisites

```bash
pip install cooler cooltools numpy scipy pandas matplotlib statsmodels bioframe
```

## Quick Start

Tell your AI agent what you want to do:
- "Compare Hi-C between treatment and control"
- "Find differential contacts"

## Example Prompts

### Basic Comparison
> "Compute log2 fold change between these two Hi-C samples"

> "Show differential contact map"

### Feature Comparison
> "Find compartment switches between conditions"

> "Compare TAD boundaries"

### Statistics
> "Test for significant differential contacts"

> "Apply FDR correction to the p-values"

## What the Agent Will Do

1. Load both cooler files
2. Normalize for sequencing depth
3. Compute log2 fold change
4. Test for significance (if replicates)
5. Visualize differential contacts

## Tips

- **Depth normalization** - Essential for fair comparison
- **Replicates** - Needed for statistical testing
- **Resolution** - Match resolution between samples
- **FDR correction** - Apply when testing many contacts


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->