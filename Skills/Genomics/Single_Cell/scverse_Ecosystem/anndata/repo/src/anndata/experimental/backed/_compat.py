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

from ..._core.xarray import Dataset2D

if TYPE_CHECKING:
    from anndata import AnnData


def has_dataset_2d(adata: AnnData) -> bool:
    if any(isinstance(annot_df, Dataset2D) for annot_df in [adata.obs, adata.var]):
        return True
    for annot_m_key in ["varm", "obsm"]:
        annot_m = getattr(adata, annot_m_key)
        if any(isinstance(maybe_df, Dataset2D) for maybe_df in annot_m.values()):
            return True
    return False

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
