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

import numpy as np
from anndata import AnnData

import scanpy as sc
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.magic]

A_list = [
    [0, 0, 7, 0, 0],
    [8, 5, 0, 2, 0],
    [6, 0, 0, 2, 5],
    [0, 0, 0, 1, 0],
    [8, 8, 2, 1, 0],
    [0, 0, 0, 4, 5],
]


def test_magic_default():
    a = np.array(A_list, dtype="float32")
    adata = AnnData(a)
    sc.external.pp.magic(adata, knn=1)
    # check raw unchanged
    np.testing.assert_array_equal(adata.raw.X, a)
    # check .X changed
    assert not np.all(a == adata.X)
    # check .X shape unchanged
    assert adata.X.shape == a.shape


def test_magic_pca_only():
    a = np.array(A_list, dtype="float32")
    # pca only
    adata = AnnData(a)
    n_pca = 3
    sc.external.pp.magic(adata, knn=1, name_list="pca_only", n_pca=n_pca)
    # check raw unchanged
    np.testing.assert_array_equal(adata.X, a)
    # check .X shape consistent with n_pca
    assert adata.obsm["X_magic"].shape == (a.shape[0], n_pca)


def test_magic_copy():
    a = np.array(A_list, dtype="float32")
    adata = AnnData(a)
    adata_copy = sc.external.pp.magic(adata, knn=1, copy=True)
    # check adata unchanged
    np.testing.assert_array_equal(adata.X, a)
    # check copy raw unchanged
    np.testing.assert_array_equal(adata_copy.raw.X, a)
    # check .X changed
    assert not np.all(a == adata_copy.X)
    # check .X shape unchanged
    assert adata_copy.X.shape == a.shape

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
