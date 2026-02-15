# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.messages import AIMessage, AIMessageChunk
from pydantic import BaseModel


class _AnyIDMixin(BaseModel):
    def __eq__(self, other: object) -> bool:
        if isinstance(other, BaseModel):
            dump = self.model_dump()
            dump.pop("id")
            other_dump = other.model_dump()
            other_dump.pop("id")
            return dump == other_dump
        return False

    __hash__ = None  # type: ignore[assignment]


class _AnyIdAIMessage(AIMessage, _AnyIDMixin):
    """AIMessage with any ID."""


class _AnyIdAIMessageChunk(AIMessageChunk, _AnyIDMixin):
    """AIMessageChunk with any ID."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
