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

# contrastiveVI

**contrastiveVI** [^ref1] (contrastive variational inference; Python class
{class}`~scvi.external.ContrastiveVI`) is a generative model for the contrastive analysis
of scRNA-seq count data that can subsequently be used for many common downstream tasks.

Contrastive analysis requires a _target_ (e.g., treated cells) and a _background_
(e.g., control cells) dataset, and contrastiveVI is designed to isolate the variations
enriched in target cells from variations shared with background cells.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/contrastiveVI_tutorial`
```

## Overview

:::{note}
This page is under construction.
:::

[^ref1]:
    Ethan Weinberger, Chris Lin, Su-In Lee (2023),
    _Isolating salient variations of interest in single-cell data with contrastiveVI_,
    [Nature Methods](https://www.nature.com/articles/s41592-023-01955-3).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->