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

# Metabolite Communication - Usage Guide

## Overview

Analyze metabolite-mediated cell-cell communication using MeboCost to infer metabolic signaling between cell types from scRNA-seq data.

## Prerequisites

```bash
pip install mebocost scanpy anndata
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze metabolite communication between my cell types"
- "Find which metabolites are secreted by macrophages"
- "Compare metabolic signaling between tumor and normal tissue"

## Example Prompts

### Basic Analysis

> "Run MeboCost on my scRNA-seq data"

> "Find significant metabolite-receptor interactions between cell types"

### Specific Metabolites

> "Which cells secrete lactate and which receive it?"

> "Analyze amino acid signaling in my dataset"

### Comparative Analysis

> "Compare metabolite communication between treatment and control"

> "Find differential metabolic signaling in tumor microenvironment"

### Visualization

> "Plot the metabolite communication network"

> "Visualize glutamine signaling flow between cell types"

## What the Agent Will Do

1. Load and verify scRNA-seq data format
2. Create MeboCost object with cell type annotations
3. Run permutation-based communication inference
4. Identify significant metabolite-receptor interactions
5. Summarize by cell type pairs and metabolites
6. Generate network visualizations

## Tips

- **Normalization** - Data must be log-normalized (scanpy standard preprocessing)
- **Cell types** - Need at least 50 cells per type for robust statistics
- **Permutations** - Use 1000 for publication; 100 for quick exploration
- **Gene symbols** - MeboCost requires gene symbols, not Ensembl IDs
- **Species** - Supports 'human' and 'mouse'


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->