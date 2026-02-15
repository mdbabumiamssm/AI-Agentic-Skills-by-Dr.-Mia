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

import scanpy as sc
import scanpy.external as sce
from testing.scanpy._helpers.data import pbmc68k_reduced
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.scanorama]


def test_scanorama_integrate():
    """Test that Scanorama integration works.

    This is a very simple test that just checks to see if the Scanorama
    integrate wrapper succesfully added a new field to ``adata.obsm``
    and makes sure it has the same dimensions as the original PCA table.
    """
    adata = pbmc68k_reduced()
    sc.pp.pca(adata)
    adata.obs["batch"] = 350 * ["a"] + 350 * ["b"]
    sce.pp.scanorama_integrate(adata, "batch", approx=False)
    assert adata.obsm["X_scanorama"].shape == adata.obsm["X_pca"].shape

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
