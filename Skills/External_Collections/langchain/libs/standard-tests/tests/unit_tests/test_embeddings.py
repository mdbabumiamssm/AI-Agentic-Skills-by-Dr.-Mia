# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any

from langchain_core.embeddings import DeterministicFakeEmbedding, Embeddings

from langchain_tests.integration_tests import EmbeddingsIntegrationTests
from langchain_tests.unit_tests import EmbeddingsUnitTests


class TestFakeEmbeddingsUnit(EmbeddingsUnitTests):
    @property
    def embeddings_class(self) -> type[Embeddings]:
        return DeterministicFakeEmbedding

    @property
    def embedding_model_params(self) -> dict[str, Any]:
        return {"size": 6}  # embedding dimension


class TestFakeEmbeddingsIntegration(EmbeddingsIntegrationTests):
    @property
    def embeddings_class(self) -> type[Embeddings]:
        return DeterministicFakeEmbedding

    @property
    def embedding_model_params(self) -> dict[str, Any]:
        return {"size": 6}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
