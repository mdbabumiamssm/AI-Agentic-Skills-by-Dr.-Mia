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

# machine-learning

## Overview

Machine learning skills for biomarker discovery, model interpretation, and predictive modeling on omics data.

**Tool type:** python | **Primary tools:** sklearn, shap, lifelines, xgboost, scvi-tools

## Skills

| Skill | Description |
|-------|-------------|
| biomarker-discovery | Boruta, mRMR, LASSO for identifying biomarkers from omics data |
| prediction-explanation | SHAP and LIME for explaining omics classifier predictions |
| model-validation | Nested CV and stratified splits for unbiased evaluation |
| atlas-mapping | scArches/scANVI for single-cell reference mapping |
| survival-analysis | Kaplan-Meier curves, Cox regression for time-to-event data |
| omics-classifiers | RandomForest, XGBoost, logistic regression for diagnostics |

## Example Prompts

- "Select biomarker genes from my expression matrix using Boruta"
- "Explain which genes are driving my classifier predictions with SHAP"
- "Run nested cross-validation on my diagnostic classifier"
- "Map my query single-cell data to the Human Lung Cell Atlas"
- "Build a survival model from my clinical data"

## Requirements

```bash
pip install scikit-learn xgboost shap lime lifelines Boruta mrmr-selection scvi-tools
```

## Related Skills

- single-cell/preprocessing - Single-cell preprocessing before ML
- differential-expression/de-results - Pre-filter genes with DE analysis
- pathway-analysis/go-enrichment - Functional enrichment of selected features


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->