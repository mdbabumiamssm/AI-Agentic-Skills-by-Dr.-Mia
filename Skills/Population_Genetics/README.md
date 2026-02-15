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

# Population Genetics Skills

Skills for population-scale genetic analysis, phylogenetics, and evolutionary biology.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **population-genetics** | Population analysis | GWAS, PCA, admixture, selection |
| **phylogenetics** | Tree analysis | ML inference, bootstrap, tree I/O |
| **comparative-genomics** | Genome comparison | Synteny, orthologs, positive selection |
| **epidemiological-genomics** | Pathogen genomics | Typing, phylodynamics, transmission |

## Key Tools

- **PLINK2** - GWAS and population analysis
- **ADMIXTURE** - Ancestry estimation
- **RAxML/IQ-TREE** - Phylogenetic inference
- **scikit-allel** - Python population genetics
- **Nextstrain** - Pathogen phylodynamics

## Example GWAS

```python
import pandas as pd
from scipy import stats

# Association testing
def run_gwas(genotypes, phenotypes):
    results = []
    for snp in genotypes.columns:
        stat, pval = stats.ttest_ind(
            phenotypes[genotypes[snp] == 0],
            phenotypes[genotypes[snp] == 1]
        )
        results.append({'snp': snp, 'pvalue': pval})
    return pd.DataFrame(results)
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->