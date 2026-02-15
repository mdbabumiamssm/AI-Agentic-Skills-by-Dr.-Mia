# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Base callback handler that can be used to handle callbacks in langchain."""

from __future__ import annotations

from langchain_core.callbacks import (
    AsyncCallbackHandler,
    BaseCallbackHandler,
    BaseCallbackManager,
    CallbackManagerMixin,
    Callbacks,
    ChainManagerMixin,
    LLMManagerMixin,
    RetrieverManagerMixin,
    RunManagerMixin,
    ToolManagerMixin,
)

__all__ = [
    "AsyncCallbackHandler",
    "BaseCallbackHandler",
    "BaseCallbackManager",
    "CallbackManagerMixin",
    "Callbacks",
    "ChainManagerMixin",
    "LLMManagerMixin",
    "RetrieverManagerMixin",
    "RunManagerMixin",
    "ToolManagerMixin",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
