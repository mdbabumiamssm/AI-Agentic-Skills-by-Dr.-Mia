# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from abc import ABC, abstractmethod

from langchain_core.runnables.config import run_in_executor
from pydantic import BaseModel, Field


class SparseVector(BaseModel, extra="forbid"):
    """Sparse vector structure."""

    indices: list[int] = Field(..., description="indices must be unique")
    values: list[float] = Field(
        ..., description="values and indices must be the same length"
    )


class SparseEmbeddings(ABC):
    """An interface for sparse embedding models to use with Qdrant."""

    @abstractmethod
    def embed_documents(self, texts: list[str]) -> list[SparseVector]:
        """Embed search docs."""

    @abstractmethod
    def embed_query(self, text: str) -> SparseVector:
        """Embed query text."""

    async def aembed_documents(self, texts: list[str]) -> list[SparseVector]:
        """Asynchronous Embed search docs."""
        return await run_in_executor(None, self.embed_documents, texts)

    async def aembed_query(self, text: str) -> SparseVector:
        """Asynchronous Embed query text."""
        return await run_in_executor(None, self.embed_query, text)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
