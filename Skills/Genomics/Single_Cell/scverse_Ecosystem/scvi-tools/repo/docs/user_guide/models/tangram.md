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

# Tangram

**Tangram** {cite:p}`Biancalani21` (Python class {class}`~scvi.external.Tangram`) maps single-cell RNA-seq data to spatial data, permitting deconvolution of cell types in spatial data like Visium.

This is a reimplementation of Tangram, which can originally be found [here](https://github.com/broadinstitute/Tangram).

## Overview

Tangram learns a matrix $M$ with shape ($n_{sc} \times n_{sp}$), in which each row sums to 1. Thus, this matrix can be viewed as a map from single cells to the spatial observations.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/tangram_scvi_tools`
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->