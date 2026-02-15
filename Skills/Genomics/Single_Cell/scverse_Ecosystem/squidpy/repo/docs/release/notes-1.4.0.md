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

# Squidpy 1.4.0 (2024-02-05)

## Bugfixes

- Fix building graph in ``knn`` and ``delaunay`` mode.
  [@michalk8](https://github.com/michalk8)
  [#792](https://github.com/scverse/squidpy/pull/792)

- Correct shuffling of annotations in ``sq.gr.nhood_enrichment``.
  [@giovp](https://github.com/giovp)
  [#775](https://github.com/scverse/squidpy/pull/775)


## Miscellaneous

- Fix napari installation.
  [@giovp](https://github.com/giovp)
  [#767](https://github.com/scverse/squidpy/pull/767)

- Made nanostring reader more flexible by adjusting loading of images.
  [@FrancescaDr](https://github.com/FrancescaDr)
  [#766](https://github.com/scverse/squidpy/pull/766)

- Fix ``sq.tl.var_by_distance`` method to support ``pandas 2.2.0``.
  [@LLehner](https://github.com/LLehner)
  [#794](https://github.com/scverse/squidpy/pull/794)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->