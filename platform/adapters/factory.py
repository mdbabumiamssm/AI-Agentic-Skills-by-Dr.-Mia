# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Dict, Any, Type
from platform.interface.llm_provider import LLMProvider
from platform.adapters.gemini_adapter import GeminiAdapter
from platform.adapters.claude_runtime_adapter import Claude37Adapter
from platform.adapters.local_adapter import LocalAdapter # Placeholder for future
import os

class LLMFactory:
    """
    Factory pattern for creating LLM Providers dynamically.
    This is the core Workflow Abstraction Layer (WAL) component.
    """
    _providers: Dict[str, Type[LLMProvider]] = {
        "gemini": GeminiAdapter,
        "anthropic": Claude37Adapter,
        "local": LocalAdapter
    }

    @staticmethod
    def create_provider(provider_type: str, config: Dict[str, Any]) -> LLMProvider:
        """
        Instantiates and initializes an LLM provider.
        """
        provider_cls = LLMFactory._providers.get(provider_type.lower())
        if not provider_cls:
            raise ValueError(f"Unknown provider type: {provider_type}")
        
        provider = provider_cls()
        success = provider.initialize(config)
        
        if not success:
            print(f"⚠️ [LLMFactory] Warning: Provider '{provider_type}' failed to initialize.")
            
        return provider

    @staticmethod
    def register_provider(name: str, provider_cls: Type[LLMProvider]):
        """
        Allows plugins to register new provider types at runtime.
        """
        LLMFactory._providers[name] = provider_cls

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
