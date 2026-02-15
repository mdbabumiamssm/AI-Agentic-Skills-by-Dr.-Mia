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

# scBasset

**scBasset** [^ref1] (Python class {class}`~scvi.external.SCBASSET`) posits a sequence-based method for representation learning of scATAC-seq data.

The advantages of ScBasset are:

-   Sequence representations allow for TF motif discovery and other sequence-based analyses.
-   scBasset is fast and scalable.

The limitations of scBasset include:

-   scBasset cannot currently leverage unobserved data and thus cannot currently be used for transfer learning tasks.

## Overview

:::{note}
This page is under construction.
:::

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/atac/scbasset`
-   {doc}`/tutorials/notebooks/atac/scbasset_batch`
```

[^ref1]:
    Yuan Han and David R. Kelley (2022),
    _scBasset: sequence-based modeling of single-cell ATAC-seq using convolutional neural networks_,
    [Nature Methods](https://www.nature.com/articles/s41592-022-01562-8).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->