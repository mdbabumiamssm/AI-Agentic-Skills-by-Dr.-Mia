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

## Preprocessing: `pp`

```{eval-rst}
.. module:: scanpy.pp
```

```{eval-rst}
.. currentmodule:: scanpy
```

Filtering of highly-variable genes, batch-effect correction, per-cell normalization, preprocessing recipes.

Any transformation of the data matrix that is not a *tool*. Other than *tools*, preprocessing steps usually don't return an easily interpretable annotation, but perform a basic transformation on the data matrix.

### Basic Preprocessing

For visual quality control, see {func}`~scanpy.pl.highest_expr_genes` and
{func}`~scanpy.pl.filter_genes_dispersion` in {mod}`scanpy.pl`.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   pp.calculate_qc_metrics
   pp.filter_cells
   pp.filter_genes
   pp.highly_variable_genes
   pp.log1p
   pp.pca
   pp.normalize_total
   pp.regress_out
   pp.scale
   pp.sample
   pp.downsample_counts
```

### Recipes

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.recipe_zheng17
   pp.recipe_weinreb17
   pp.recipe_seurat
```

### Batch effect correction

Also see {ref}`data-integration`. Note that a simple batch correction method is available via {func}`pp.regress_out`. Checkout {mod}`scanpy.external` for more.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.combat
```

### Doublet detection

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.scrublet
   pp.scrublet_simulate_doublets
```

### Neighbors

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.neighbors

```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->