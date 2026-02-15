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
from typing_extensions import override


class SequentialRetriever(BaseRetriever):
    """Test util that returns a sequence of documents."""

    sequential_responses: list[list[Document]]
    response_index: int = 0

    @override
    def _get_relevant_documents(
        self,
        query: str,
        **kwargs: Any,
    ) -> list[Document]:
        if self.response_index >= len(self.sequential_responses):
            return []
        self.response_index += 1
        return self.sequential_responses[self.response_index - 1]

    @override
    async def _aget_relevant_documents(
        self,
        query: str,
        **kwargs: Any,
    ) -> list[Document]:
        return self._get_relevant_documents(query)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
