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

import pytest

import scanpy as sc


# TODO: reset everything
@pytest.mark.parametrize("_attempt", range(3))
def test_resets(_attempt: int) -> None:
    """Test that changes made reset."""
    assert sc.settings.autoshow
    sc.settings.autoshow = False
    assert not sc.settings.autoshow


def test_set_figure_params_warns() -> None:
    with pytest.warns(FutureWarning, match=r"scanpy\.set_figure_params"):
        sc.settings.set_figure_params()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
