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
from testing.scanpy._helpers.data import pbmc3k
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.wishbone]


def test_run_wishbone():
    adata = pbmc3k()
    sc.pp.normalize_per_cell(adata)
    sc.pp.neighbors(adata, n_pcs=15, n_neighbors=10)
    sc.pp.pca(adata)
    sc.tl.tsne(adata=adata, n_pcs=5, perplexity=30)
    sc.tl.diffmap(adata, n_comps=10)

    sce.tl.wishbone(
        adata=adata,
        start_cell="ACAAGAGACTTATC-1",
        components=[2, 3],
        num_waypoints=150,
    )
    assert all(k in adata.obs for k in ["trajectory_wishbone", "branch_wishbone"]), (
        "Run Wishbone Error!"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
