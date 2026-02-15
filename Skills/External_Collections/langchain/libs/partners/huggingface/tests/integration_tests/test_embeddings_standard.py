# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test HuggingFace embeddings."""

from langchain_tests.integration_tests import EmbeddingsIntegrationTests

from langchain_huggingface.embeddings import (
    HuggingFaceEmbeddings,
    HuggingFaceEndpointEmbeddings,
)


class TestHuggingFaceEmbeddings(EmbeddingsIntegrationTests):
    @property
    def embeddings_class(self) -> type[HuggingFaceEmbeddings]:
        return HuggingFaceEmbeddings

    @property
    def embedding_model_params(self) -> dict:
        return {"model_name": "sentence-transformers/all-mpnet-base-v2"}


class TestHuggingFaceEndpointEmbeddings(EmbeddingsIntegrationTests):
    @property
    def embeddings_class(self) -> type[HuggingFaceEndpointEmbeddings]:
        return HuggingFaceEndpointEmbeddings

    @property
    def embedding_model_params(self) -> dict:
        return {"model": "sentence-transformers/all-mpnet-base-v2"}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
