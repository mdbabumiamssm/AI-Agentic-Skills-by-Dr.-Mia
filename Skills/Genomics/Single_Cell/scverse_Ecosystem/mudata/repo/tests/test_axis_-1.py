# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import unittest

import numpy as np
import pytest
from anndata import AnnData

from mudata import MuData


@pytest.mark.usefixtures("filepath_h5mu")
class TestMuData:
    def test_create(self):
        n, d_raw, d_preproc = 100, 900, 300

        a_raw = AnnData(np.random.normal(size=(n, d_raw)))
        a_preproc = a_raw[
            :, np.sort(np.random.choice(np.arange(d_raw), d_preproc, replace=False))
        ].copy()

        mdata = MuData({"raw": a_raw, "preproc": a_preproc}, axis=-1)

        assert mdata.n_obs == n
        assert mdata.n_vars == d_raw


if __name__ == "__main__":
    unittest.main()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
