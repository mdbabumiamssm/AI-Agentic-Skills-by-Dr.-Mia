# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import numpy as np
import pytest

from scvi.data import synthetic_iid


@pytest.mark.parametrize("sparse_format", ["csr_matrix", "csc_matrix", None])
def test_synthetic_iid_sparse_format(sparse_format: str | None):
    _ = synthetic_iid(sparse_format=sparse_format)


@pytest.mark.parametrize("return_mudata", [False, True])
def test_synthetic_iid_coords(return_mudata: bool):
    adata = synthetic_iid(return_mudata=return_mudata, generate_coordinates=True)
    assert "coordinates" in adata.obsm
    assert isinstance(adata.obsm["coordinates"], np.ndarray)
    assert adata.obsm["coordinates"].shape == (adata.n_obs, 2)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
