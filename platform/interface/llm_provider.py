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
from typing import List, Dict, Any, Optional
from platform.schema.io_types import LLMRequest, LLMResponse

class LLMProvider(ABC):
    """
    Abstract Base Class for all LLM providers (Gemini, OpenAI, Anthropic, Local).
    Enforces a standardized interface for the BioKernel.
    """

    @abstractmethod
    def initialize(self, config: Dict[str, Any]) -> bool:
        """
        Initialize the provider with configuration (API keys, model names, etc.).
        Returns True if successful.
        """
        pass

    @abstractmethod
    async def generate(self, request: LLMRequest) -> LLMResponse:
        """
        Standard generation method. 
        """
        pass

    @abstractmethod
    async def run_reasoning_loop(self, goal: str, tools: List[Dict]) -> str:
        """
        Executes a ReAct or Chain-of-Thought loop.
        """
        pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
