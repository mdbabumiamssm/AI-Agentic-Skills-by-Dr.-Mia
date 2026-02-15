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

# Polygenic Risk Scores - Usage Guide

## Overview

Calculate polygenic risk scores (PRS) from GWAS summary statistics using PRSice-2 (clumping/thresholding), LDpred2 (Bayesian shrinkage), or PRS-CS to predict disease risk.

## Prerequisites

```bash
# PRSice-2
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip

# PRS-CS
git clone https://github.com/getian107/PRScs.git
pip install scipy h5py
```

```r
# LDpred2
install.packages('bigsnpr')
```

## Quick Start

Tell your AI agent:
- "Calculate PRS for coronary artery disease using this GWAS"
- "Run LDpred2-auto on my genotype data"
- "Normalize my PRS scores to Z-scores"
- "Stratify samples into risk categories"

## Example Prompts

### Basic PRS

> "Run PRSice-2 on my GWAS summary stats and genotype data"

> "Calculate PRS at multiple p-value thresholds"

### Advanced Methods

> "Use LDpred2-auto for PRS calculation"

> "Run PRS-CS with the European LD reference panel"

### Interpretation

> "Convert my raw PRS to population percentiles"

> "Categorize samples into low/average/high risk groups"

## What the Agent Will Do

1. Match GWAS variants to target genotype data
2. Apply clumping to remove LD-correlated variants (PRSice-2)
3. Or apply Bayesian shrinkage (LDpred2/PRS-CS)
4. Calculate weighted sum of effect alleles
5. Normalize scores and compute percentiles
6. Validate using phenotype if available

## Method Comparison

| Method | Approach | Best For |
|--------|----------|----------|
| PRSice-2 | Clump + threshold | Quick analysis, multiple thresholds |
| LDpred2-auto | Bayesian, auto-tuned | Best accuracy, research |
| PRS-CS | Bayesian continuous shrinkage | External LD panel |

## Tips

- **GWAS and target must match** ancestry or adjust with PCA
- **LDpred2-auto** is recommended for best accuracy without tuning
- **Use external LD panel** when target sample is small
- **Normalize PRS** to Z-scores for interpretability
- **PGS Catalog** provides pre-computed weights for many traits
- **Validate with phenotype** when available (AUC, R2)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->