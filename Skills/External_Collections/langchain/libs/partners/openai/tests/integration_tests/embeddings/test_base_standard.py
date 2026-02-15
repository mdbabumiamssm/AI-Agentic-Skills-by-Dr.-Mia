# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Standard LangChain interface tests"""

from langchain_core.embeddings import Embeddings
from langchain_tests.integration_tests.embeddings import EmbeddingsIntegrationTests

from langchain_openai import OpenAIEmbeddings


class TestOpenAIStandard(EmbeddingsIntegrationTests):
    @property
    def embeddings_class(self) -> type[Embeddings]:
        return OpenAIEmbeddings

    @property
    def embedding_model_params(self) -> dict:
        return {"model": "text-embedding-3-small", "dimensions": 128}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
