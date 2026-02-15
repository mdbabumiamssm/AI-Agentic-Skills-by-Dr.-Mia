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
import pytest

import scanpy as sc


def test_sim_toggleswitch() -> None:
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata_sim = sc.tl.sim("toggleswitch")
    with pytest.warns(UserWarning, match=r"Observation names are not unique"):
        adata_ds = sc.datasets.toggleswitch()
    np.allclose(adata_sim.X, adata_ds.X, np.finfo(np.float32).eps)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
