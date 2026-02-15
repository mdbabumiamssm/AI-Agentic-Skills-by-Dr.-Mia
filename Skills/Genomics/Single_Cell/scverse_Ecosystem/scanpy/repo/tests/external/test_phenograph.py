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
import pandas as pd
from anndata import AnnData

import scanpy as sc
import scanpy.external as sce
from testing.scanpy._pytest.marks import needs

pytestmark = [needs.phenograph]


def test_phenograph():
    df = np.random.rand(1000, 40)
    dframe = pd.DataFrame(df)
    dframe.index, dframe.columns = (map(str, dframe.index), map(str, dframe.columns))
    adata = AnnData(dframe)
    sc.pp.pca(adata, n_comps=20)
    sce.tl.phenograph(adata, clustering_algo="leiden", k=50)
    assert adata.obs["pheno_leiden"].shape[0], "phenograph_Community Detection Error!"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
