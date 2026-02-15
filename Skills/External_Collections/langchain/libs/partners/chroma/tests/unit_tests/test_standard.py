# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from collections.abc import Generator

import pytest
from langchain_core.vectorstores import VectorStore
from langchain_tests.integration_tests.vectorstores import VectorStoreIntegrationTests
from pytest_benchmark.fixture import BenchmarkFixture  # type: ignore[import-untyped]

from langchain_chroma import Chroma


class TestChromaStandard(VectorStoreIntegrationTests):
    @pytest.fixture
    def vectorstore(self) -> Generator[VectorStore, None, None]:  # type: ignore[override]
        """Get an empty vectorstore for unit tests."""
        store = Chroma(embedding_function=self.get_embeddings())
        try:
            yield store
        finally:
            store.delete_collection()


@pytest.mark.benchmark
def test_chroma_init_time(benchmark: BenchmarkFixture) -> None:
    """Test Chroma initialization time."""

    def _init_chroma() -> None:
        for _ in range(10):
            Chroma()

    benchmark(_init_chroma)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
