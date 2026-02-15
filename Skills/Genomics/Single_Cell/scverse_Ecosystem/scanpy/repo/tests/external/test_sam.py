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

import scanpy as sc
import scanpy.external as sce
from testing.scanpy._helpers.data import pbmc3k
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.samalg]


def test_sam():
    adata_ref = pbmc3k()
    ix = np.random.choice(adata_ref.shape[0], size=200, replace=False)
    adata = adata_ref[ix, :].copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sce.tl.sam(adata, inplace=True)
    uns_keys = list(adata.uns.keys())
    obsm_keys = list(adata.obsm.keys())
    assert all(["sam" in uns_keys, "X_umap" in obsm_keys, "neighbors" in uns_keys])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
