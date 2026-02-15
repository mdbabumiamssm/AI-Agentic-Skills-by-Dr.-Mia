# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Standard unit tests for ExaSearchRetriever."""

import pytest
from pytest_benchmark.fixture import BenchmarkFixture  # type: ignore[import-untyped]

from langchain_exa import ExaSearchRetriever


@pytest.mark.benchmark
def test_exa_retriever_init_time(benchmark: BenchmarkFixture) -> None:
    """Test ExaSearchRetriever initialization time."""

    def _init_exa_retriever() -> None:
        for _ in range(10):
            ExaSearchRetriever()

    benchmark(_init_exa_retriever)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
