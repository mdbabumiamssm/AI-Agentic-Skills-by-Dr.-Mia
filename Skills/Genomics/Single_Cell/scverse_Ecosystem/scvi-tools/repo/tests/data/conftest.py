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


@pytest.fixture
def adata():
    anndata = synthetic_iid()
    raw_counts = anndata.X.copy()
    anndata.layers["raw"] = raw_counts
    anndata.obs["cont1"] = np.random.normal(size=(anndata.shape[0],))
    anndata.obs["cont2"] = np.random.normal(size=(anndata.shape[0],))
    anndata.obs["cat1"] = np.random.randint(0, 5, size=(anndata.shape[0],))
    anndata.obs["cat2"] = np.random.randint(0, 5, size=(anndata.shape[0],))

    return anndata


adata1 = adata2 = adata

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
