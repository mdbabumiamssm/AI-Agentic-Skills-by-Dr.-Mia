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

from typing import Annotated, Any, Literal

from pydantic import Field

from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings


class GoogleAIPromptExecutionSettings(PromptExecutionSettings):
    """Google AI Prompt Execution Settings."""

    stop_sequences: Annotated[list[str] | None, Field(max_length=5)] = None
    response_mime_type: Literal["text/plain", "application/json"] | None = None
    response_schema: Any | None = None
    candidate_count: Annotated[int | None, Field(ge=1)] = None
    max_output_tokens: Annotated[int | None, Field(ge=1)] = None
    temperature: Annotated[float | None, Field(ge=0.0, le=2.0)] = None
    top_p: float | None = None
    top_k: int | None = None


class GoogleAITextPromptExecutionSettings(GoogleAIPromptExecutionSettings):
    """Google AI Text Prompt Execution Settings."""

    pass


class GoogleAIChatPromptExecutionSettings(GoogleAIPromptExecutionSettings):
    """Google AI Chat Prompt Execution Settings."""

    tools: Annotated[
        list[dict[str, Any]] | None,
        Field(
            description="Do not set this manually. It is set by the service based "
            "on the function choice configuration.",
        ),
    ] = None
    tool_config: Annotated[
        dict[str, Any] | None,
        Field(
            description="Do not set this manually. It is set by the service based "
            "on the function choice configuration.",
        ),
    ] = None


class GoogleAIEmbeddingPromptExecutionSettings(PromptExecutionSettings):
    """Google AI Embedding Prompt Execution Settings."""

    output_dimensionality: Annotated[int | None, Field(le=768)] = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
