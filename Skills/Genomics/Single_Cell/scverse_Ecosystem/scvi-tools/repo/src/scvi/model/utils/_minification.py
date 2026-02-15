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

from scipy.sparse import csr_matrix

from scvi import REGISTRY_KEYS

if TYPE_CHECKING:
    from anndata import AnnData
    from mudata import MuData

    from scvi.data import AnnDataManager


def get_minified_adata_scrna(
    adata_manager: AnnDataManager,
    keep_count_data: bool = False,
) -> AnnData:
    """Returns a minified AnnData.

    Parameters
    ----------
    adata_manager
        Manager with original AnnData, of which we want to create a minified version.
    keep_count_data
        If True, the count data is kept in the minified data. If False, the count data is removed.
    """
    adata = adata_manager.adata.copy()
    if keep_count_data:
        pass
    else:
        del adata.raw
        counts = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        all_zeros = csr_matrix(counts.shape)
        X = all_zeros
        layers = {layer: all_zeros.copy() for layer in adata_manager.adata.layers}
        adata.X = X
        adata.layers = layers
    return adata


def get_minified_mudata(
    adata_manager: AnnDataManager,
    keep_count_data: bool = False,
) -> MuData:
    """Returns a minified MuData that works for most multi modality models (MULTIVI, TOTALVI).

    Parameters
    ----------
    adata_manager
        Manager with original MuData, of which we want to create a minified version.
    keep_count_data
        If True, the count data is kept in the minified data. If False, the count data is removed.
    """
    mdata = adata_manager.adata.copy()
    if keep_count_data:
        pass
    else:
        for modality in mdata.mod_names:
            del mdata[modality].raw
            all_zeros = csr_matrix(mdata[modality].X.shape)
            mdata[modality].X = all_zeros.copy()
            if len(mdata[modality].layers) > 0:
                layers = {layer: all_zeros.copy() for layer in mdata[modality].layers}
                mdata[modality].layers = layers
    return mdata

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
