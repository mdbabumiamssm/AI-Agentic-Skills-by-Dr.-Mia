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
from langchain_core.vectorstores import (
    InMemoryVectorStore,
    VectorStore,
)

from langchain_tests.integration_tests.vectorstores import VectorStoreIntegrationTests


class TestInMemoryVectorStore(VectorStoreIntegrationTests):
    @pytest.fixture
    def vectorstore(self) -> VectorStore:
        embeddings = self.get_embeddings()
        return InMemoryVectorStore(embedding=embeddings)


class WithoutGetByIdsVectorStore(InMemoryVectorStore):
    """InMemoryVectorStore that does not implement get_by_ids."""

    get_by_ids = VectorStore.get_by_ids


class TestWithoutGetByIdVectorStore(VectorStoreIntegrationTests):
    @pytest.fixture
    def vectorstore(self) -> VectorStore:
        embeddings = self.get_embeddings()
        return WithoutGetByIdsVectorStore(embedding=embeddings)

    @property
    def has_get_by_ids(self) -> bool:
        return False

    def test_get_by_ids_fails(self, vectorstore: VectorStore) -> None:
        with pytest.raises(
            NotImplementedError,
            match="WithoutGetByIdsVectorStore does not yet support get_by_ids",
        ):
            vectorstore.get_by_ids(["id1", "id2"])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
