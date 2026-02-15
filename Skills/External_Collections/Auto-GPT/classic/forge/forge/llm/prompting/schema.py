# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import enum

from pydantic import BaseModel, Field

from forge.llm.providers.schema import (
    ChatMessage,
    ChatMessageDict,
    CompletionModelFunction,
)


class LanguageModelClassification(str, enum.Enum):
    """The LanguageModelClassification is a functional description of the model.

    This is used to determine what kind of model to use for a given prompt.
    Sometimes we prefer a faster or cheaper model to accomplish a task when
    possible.
    """

    FAST_MODEL = "fast_model"
    SMART_MODEL = "smart_model"


class ChatPrompt(BaseModel):
    messages: list[ChatMessage]
    functions: list[CompletionModelFunction] = Field(default_factory=list)
    prefill_response: str = ""

    def raw(self) -> list[ChatMessageDict]:
        return [m.model_dump() for m in self.messages]  # type: ignore

    def __str__(self):
        return "\n\n".join(
            f"{m.role.value.upper()}: {m.content}" for m in self.messages
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
