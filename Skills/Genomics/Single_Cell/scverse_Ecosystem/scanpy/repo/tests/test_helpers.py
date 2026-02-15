# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test our test helpers."""

from __future__ import annotations

import numpy as np

from testing.scanpy._helpers import random_mask


def test_random_mask():
    ns_true = np.array([int(random_mask(4).sum()) for _ in range(1000)])
    np.testing.assert_equal(ns_true, [2] * 1000)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
