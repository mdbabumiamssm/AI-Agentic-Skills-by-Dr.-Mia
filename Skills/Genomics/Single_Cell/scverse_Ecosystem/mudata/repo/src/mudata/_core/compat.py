# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping

    from anndata import AnnData
    from anndata._core.raw import Raw

try:
    from anndata._core.aligned_mapping import AlignedView, AxisArrays, PairwiseArrays
except ImportError:
    # anndata < 0.10.9
    from anndata._core.aligned_mapping import (
        AlignedViewMixin as AlignedView,
    )
    from anndata._core.aligned_mapping import (
        AxisArrays as AxisArraysLegacy,
    )
    from anndata._core.aligned_mapping import (
        AxisArraysBase,
    )
    from anndata._core.aligned_mapping import (
        PairwiseArrays as PairwiseArraysLegacy,
    )

    class AxisArrays(AxisArraysLegacy):
        def __init__(
            self,
            parent: AnnData | Raw,
            axis: int,
            store: Mapping | AxisArraysBase | None = None,
        ):
            super().__init__(parent, axis=axis, vals=store)

    class PairwiseArrays(PairwiseArraysLegacy):
        def __init__(
            self,
            parent: AnnData,
            axis: int,
            store: Mapping | None = None,
        ):
            super().__init__(parent, axis=axis, vals=store)


__all__ = ["AlignedView", "AxisArrays", "PairwiseArrays"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
