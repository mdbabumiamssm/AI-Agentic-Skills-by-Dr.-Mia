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

# Experimental Design Skills

Skills for planning and optimizing biological experiments.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **experimental-design** | Study design | Power analysis, sample size, batch effects |

## Key Capabilities

- **Power Analysis** - Calculate required sample sizes
- **Sample Size Estimation** - For various study types
- **Batch Design** - Minimize batch effects
- **Randomization** - Treatment allocation
- **Multiple Testing** - Correction methods

## Key Tools

- **statsmodels** - Statistical power calculations
- **scipy.stats** - Statistical distributions
- **pingouin** - Statistical analysis
- **sva** - Batch effect correction (R)

## Example Power Analysis

```python
from statsmodels.stats.power import TTestIndPower

# Power analysis for t-test
analysis = TTestIndPower()

# Calculate required sample size
sample_size = analysis.solve_power(
    effect_size=0.5,  # Cohen's d
    power=0.8,
    alpha=0.05,
    ratio=1.0,  # equal group sizes
    alternative='two-sided'
)

print(f"Required samples per group: {sample_size:.0f}")
```

## Batch Design Strategy

```python
def balanced_batch_design(samples, n_batches, covariates):
    """
    Create balanced batch assignments minimizing confounding.
    """
    from sklearn.model_selection import StratifiedKFold

    # Stratify by key covariates
    strata = covariates['condition'].astype(str) + '_' + \
             covariates['sex'].astype(str)

    skf = StratifiedKFold(n_splits=n_batches, shuffle=True)
    batches = np.zeros(len(samples))

    for batch_idx, (_, indices) in enumerate(skf.split(samples, strata)):
        batches[indices] = batch_idx

    return batches
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->