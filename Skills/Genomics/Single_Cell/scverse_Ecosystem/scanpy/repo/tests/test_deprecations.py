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
from testing.scanpy._helpers.data import pbmc68k_reduced


def test_deprecate_multicore_tsne() -> None:
    pbmc = pbmc68k_reduced()

    with pytest.warns(
        UserWarning, match=r"calling tsne with `n_jobs` > 1 would use MulticoreTSNE"
    ):
        sc.tl.tsne(pbmc, n_jobs=2)

    with (
        pytest.warns(FutureWarning, match=r"Argument `use_fast_tsne` is deprecated"),
        pytest.warns(ImportWarning, match=r"MulticoreTSNE"),
    ):
        sc.tl.tsne(pbmc, use_fast_tsne=True)


def test_deprecate_use_highly_variable_genes():
    pbmc = pbmc68k_reduced()

    with pytest.warns(
        FutureWarning, match="Argument `use_highly_variable` is deprecated"
    ):
        sc.pp.pca(pbmc, use_highly_variable=True)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
