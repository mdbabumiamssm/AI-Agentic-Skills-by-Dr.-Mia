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


class OllamaPromptExecutionSettings(PromptExecutionSettings):
    """Settings for Ollama prompt execution."""

    format: Literal["json"] | dict[str, Any] | None = None
    options: dict[str, Any] | None = None


class OllamaTextPromptExecutionSettings(OllamaPromptExecutionSettings):
    """Settings for Ollama text prompt execution."""

    system: str | None = None
    template: str | None = None
    context: str | None = None
    raw: bool | None = None


class OllamaChatPromptExecutionSettings(OllamaPromptExecutionSettings):
    """Settings for Ollama chat prompt execution."""

    tools: Annotated[
        list[dict[str, Any]] | None,
        Field(
            description="Do not set this manually. It is set by the service based "
            "on the function choice configuration.",
        ),
    ] = None


class OllamaEmbeddingPromptExecutionSettings(OllamaPromptExecutionSettings):
    """Settings for Ollama embedding prompt execution."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
