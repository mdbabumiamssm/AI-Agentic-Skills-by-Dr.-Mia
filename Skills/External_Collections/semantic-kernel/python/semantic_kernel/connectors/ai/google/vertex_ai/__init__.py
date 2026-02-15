# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

import warnings

from semantic_kernel.connectors.ai.google.vertex_ai.services.vertex_ai_chat_completion import VertexAIChatCompletion
from semantic_kernel.connectors.ai.google.vertex_ai.services.vertex_ai_text_completion import VertexAITextCompletion
from semantic_kernel.connectors.ai.google.vertex_ai.services.vertex_ai_text_embedding import VertexAITextEmbedding
from semantic_kernel.connectors.ai.google.vertex_ai.vertex_ai_prompt_execution_settings import (
    VertexAIChatPromptExecutionSettings,
    VertexAIEmbeddingPromptExecutionSettings,
    VertexAIPromptExecutionSettings,
    VertexAITextPromptExecutionSettings,
)

# Deprecation warning for the entire Vertex AI package
warnings.warn(
    "The `semantic_kernel.connectors.ai.google.vertex_ai` package is deprecated and will be removed after 01/01/2026. "
    "Please use `semantic_kernel.connectors.ai.google` instead for Google AI services.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "VertexAIChatCompletion",
    "VertexAIChatPromptExecutionSettings",
    "VertexAIEmbeddingPromptExecutionSettings",
    "VertexAIPromptExecutionSettings",
    "VertexAITextCompletion",
    "VertexAITextEmbedding",
    "VertexAITextPromptExecutionSettings",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
