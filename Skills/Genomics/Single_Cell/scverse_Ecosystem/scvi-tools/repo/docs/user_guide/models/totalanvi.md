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

# TotalANVI

:::{note}
This page is under construction.
:::

**TotalANVI** [^ref1] (Python class {class}`~scvi.external.TOTALANVI`) is a semi-supervised generative model of CITE-seq RNA and protein data.
Similar to how scANVI extends scVI, TotalANVI can be treated as an extension of TotalVI that can leverage cell type annotations
for a subset of the cells present in the data sets to infer the states of the rest of the cells as well as impute missing proteins expression

The advantages of TotalANVI are:

-   Comprehensive in capabilities.
-   Scalable to very large datasets (>1 million cells).

The limitations of TotalANVI include:

-   Effectively requires a GPU for fast inference.
-   May not scale to a very large number of cell types.

```{topic} Tutorials:

-   Work in progress.
```

## Preliminaries

## Generative process

## Inference

## Training details

## Tasks

### Cell type label prediction


[^ref1]:


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->