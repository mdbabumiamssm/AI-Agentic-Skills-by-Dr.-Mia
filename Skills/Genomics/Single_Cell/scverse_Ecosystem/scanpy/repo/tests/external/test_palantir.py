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

import scanpy.external as sce
from testing.scanpy._helpers.data import pbmc3k_processed
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.palantir]


def test_palantir_core():
    adata = pbmc3k_processed()

    sce.tl.palantir(adata=adata, n_components=5, knn=30)
    assert adata.layers["palantir_imp"].shape[0], "palantir_imp matrix Error!"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
