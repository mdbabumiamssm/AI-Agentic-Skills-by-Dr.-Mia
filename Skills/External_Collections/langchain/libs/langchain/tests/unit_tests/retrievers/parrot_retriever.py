# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.documents import Document
from langchain_core.retrievers import BaseRetriever


class FakeParrotRetriever(BaseRetriever):
    """Test util that parrots the query back as documents."""

    def _get_relevant_documents(  # type: ignore[override]
        self,
        query: str,
    ) -> list[Document]:
        return [Document(page_content=query)]

    async def _aget_relevant_documents(  # type: ignore[override]
        self,
        query: str,
    ) -> list[Document]:
        return [Document(page_content=query)]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
