# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.embeddings import Embeddings
from langchain_tests.unit_tests.embeddings import EmbeddingsUnitTests

from langchain_openai import AzureOpenAIEmbeddings


class TestAzureOpenAIStandard(EmbeddingsUnitTests):
    @property
    def embeddings_class(self) -> type[Embeddings]:
        return AzureOpenAIEmbeddings

    @property
    def embedding_model_params(self) -> dict:
        return {"api_key": "api_key", "azure_endpoint": "https://endpoint.com"}

    @property
    def init_from_env_params(self) -> tuple[dict, dict, dict]:
        return (
            {
                "AZURE_OPENAI_API_KEY": "api_key",
                "AZURE_OPENAI_ENDPOINT": "https://endpoint.com",
                "AZURE_OPENAI_AD_TOKEN": "token",
                "OPENAI_ORG_ID": "org_id",
                "OPENAI_API_VERSION": "yyyy-mm-dd",
                "OPENAI_API_TYPE": "type",
            },
            {},
            {
                "openai_api_key": "api_key",
                "azure_endpoint": "https://endpoint.com",
                "azure_ad_token": "token",
                "openai_organization": "org_id",
                "openai_api_version": "yyyy-mm-dd",
                "openai_api_type": "type",
            },
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
