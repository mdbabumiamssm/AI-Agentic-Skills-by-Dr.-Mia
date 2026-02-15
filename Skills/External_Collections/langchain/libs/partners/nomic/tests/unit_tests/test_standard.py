# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Unit tests for standard tests in Nomic partner integration."""

import pytest
from pytest_benchmark.fixture import BenchmarkFixture  # type: ignore[import]

from langchain_nomic import NomicEmbeddings


@pytest.mark.benchmark
def test_nomic_embeddings_init_time(benchmark: BenchmarkFixture) -> None:
    """Test NomicEmbeddings initialization time."""

    def _init_nomic_embeddings() -> None:
        for _ in range(10):
            NomicEmbeddings(model="test")

    benchmark(_init_nomic_embeddings)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
