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

# Squidpy 1.1.1 (2021-08-16)

## Features

- Allow defining cylindrical shells in {func}`squidpy.gr.spatial_neighbors` by using the ``radius`` argument.
  Also rename ``n_neigh``, ``n_neigh_grid`` arguments to ``n_neighs``.
  [#393](https://github.com/scverse/squidpy/pull/393)

- Allow specifying gene symbols from :attr:`anndata.AnnData.var` in {func}`squidpy.gr.ligrec`.
  [#395](https://github.com/scverse/squidpy/pull/395)


## Bugfixes

- Fix sometimes incorrectly transposing dimensions when reading TIFF files.
  [#390](https://github.com/scverse/squidpy/pull/390)


## Miscellaneous

- Increase performance of Delaunay graph creation in {func}`squidpy.gr.spatial_neighbors`.
  [#381](https://github.com/scverse/squidpy/pull/381)

- Update ``mypy`` type-checking and use PEP 604 for type annotations.
  [#396](https://github.com/scverse/squidpy/pull/396)


## Documentation

- Enable ``towncrier`` for release notes generation.
  [#397](https://github.com/scverse/squidpy/pull/397)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->