# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test Ollama embeddings."""

import os

from langchain_tests.integration_tests import EmbeddingsIntegrationTests

from langchain_ollama.embeddings import OllamaEmbeddings

MODEL_NAME = os.environ.get("OLLAMA_TEST_MODEL", "llama3.1")


class TestOllamaEmbeddings(EmbeddingsIntegrationTests):
    @property
    def embeddings_class(self) -> type[OllamaEmbeddings]:
        return OllamaEmbeddings

    @property
    def embedding_model_params(self) -> dict:
        return {"model": MODEL_NAME}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
