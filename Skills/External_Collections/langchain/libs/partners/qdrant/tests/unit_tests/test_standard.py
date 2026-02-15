# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest
from langchain_core.embeddings import Embeddings
from pytest_benchmark.fixture import BenchmarkFixture  # type: ignore[import-untyped]

from langchain_qdrant import QdrantVectorStore


class MockEmbeddings(Embeddings):
    """Mock embeddings for testing."""

    def embed_documents(self, texts: list[str]) -> list[list[float]]:
        """Mock embed_documents method."""
        return [[1.0, 2.0, 3.0] for _ in texts]

    def embed_query(self) -> list[float]:  # type: ignore[override]
        """Mock embed_query method."""
        return [1.0, 2.0, 3.0]


@pytest.mark.benchmark
def test_qdrant_vectorstore_init_time(benchmark: BenchmarkFixture) -> None:
    """Test QdrantVectorStore initialization time."""

    def _init_qdrant_vectorstore() -> None:
        for _ in range(10):
            QdrantVectorStore.from_texts(
                texts=["test"],
                embedding=MockEmbeddings(),
                location=":memory:",
                collection_name="test",
            )

    benchmark(_init_qdrant_vectorstore)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
