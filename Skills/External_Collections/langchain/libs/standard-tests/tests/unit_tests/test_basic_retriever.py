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

from langchain_core.documents import Document
from langchain_core.retrievers import BaseRetriever

from langchain_tests.integration_tests import RetrieversIntegrationTests


class ParrotRetriever(BaseRetriever):
    parrot_name: str
    k: int = 3

    def _get_relevant_documents(self, query: str, **kwargs: Any) -> list[Document]:
        k = kwargs.get("k", self.k)
        return [Document(page_content=f"{self.parrot_name} says: {query}")] * k


class TestParrotRetrieverIntegration(RetrieversIntegrationTests):
    @property
    def retriever_constructor(self) -> type[ParrotRetriever]:
        return ParrotRetriever

    @property
    def retriever_constructor_params(self) -> dict[str, Any]:
        return {"parrot_name": "Polly"}

    @property
    def retriever_query_example(self) -> str:
        return "parrot"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
