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

# Off-Target Prediction - Usage Guide

## Overview

Predict and score potential off-target cleavage sites for CRISPR guides using Cas-OFFinder and CFD scoring algorithms.

## Prerequisites

```bash
# Cas-OFFinder (download binary from http://www.rgenome.net/cas-offinder/)
# Or compile from source:
git clone https://github.com/nickloman/cas-offinder
cd cas-offinder && cmake . && make

pip install pandas biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Check off-target sites for my guide sequence ATCGATCGATCGATCGATCG"
- "Find potential off-targets for these 5 guides in the human genome"
- "Calculate specificity scores for my CRISPR guides"

## Example Prompts

### Single Guide Analysis

> "Find off-targets for ATCGATCGATCGATCGATCG allowing up to 4 mismatches"

> "What is the specificity score for this guide sequence?"

### Batch Analysis

> "Rank my 10 guides by off-target specificity"

> "Find the most specific guide from this list"

### Context-Aware

> "Check if any off-targets fall in coding regions"

> "Flag off-targets in tumor suppressor genes"

## What the Agent Will Do

1. Accept guide sequence(s) and genome reference
2. Run Cas-OFFinder to identify potential off-target sites
3. Calculate CFD scores for each off-target
4. Compute aggregate specificity scores
5. Flag high-risk off-targets (low mismatch count, high CFD)
6. Return ranked results with recommendations

## Tips

- **Mismatch tolerance** - 4 mismatches balances speed and sensitivity; use 3 for faster results
- **CFD threshold** - Sites with CFD > 0.5 have significant cleavage risk
- **Exonic off-targets** - Prioritize avoiding off-targets in coding regions
- **GPU acceleration** - Use GPU mode for faster genome-wide searches
- **Multiple guides** - Always compare specificity across candidate guides


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->