# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import abc
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from forge.llm.providers import AssistantChatMessage

from .schema import ChatPrompt, LanguageModelClassification


class PromptStrategy(abc.ABC):
    @property
    @abc.abstractmethod
    def llm_classification(self) -> LanguageModelClassification:
        ...

    @abc.abstractmethod
    def build_prompt(self, *_, **kwargs) -> ChatPrompt:
        ...

    @abc.abstractmethod
    def parse_response_content(self, response: "AssistantChatMessage") -> Any:
        ...

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
