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
from langchain_tests.unit_tests.embeddings import EmbeddingsUnitTests

from langchain_openai import OpenAIEmbeddings


class TestOpenAIStandard(EmbeddingsUnitTests):
    @property
    def embeddings_class(self) -> type[Embeddings]:
        return OpenAIEmbeddings

    @property
    def init_from_env_params(self) -> tuple[dict, dict, dict]:
        return (
            {
                "OPENAI_API_KEY": "api_key",
                "OPENAI_ORG_ID": "org_id",
                "OPENAI_API_BASE": "api_base",
                "OPENAI_PROXY": "https://proxy.com",
            },
            {},
            {
                "openai_api_key": "api_key",
                "openai_organization": "org_id",
                "openai_api_base": "api_base",
                "openai_proxy": "https://proxy.com",
            },
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
