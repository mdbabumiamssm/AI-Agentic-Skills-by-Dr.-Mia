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

import logging
from typing import Annotated, Literal

from pydantic import Field

from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings

logger = logging.getLogger(__name__)


class OpenAITextToAudioExecutionSettings(PromptExecutionSettings):
    """Request settings for OpenAI text to audio services."""

    ai_model_id: str | None = Field(None, serialization_alias="model")
    input: str | None = Field(
        None, description="Do not set this manually. It is set by the service based on the text content."
    )
    voice: Literal["alloy", "ash", "ballad", "coral", "echo", "fable", "onyx", "nova", "sage", "shimmer"] = "alloy"
    response_format: Literal["mp3", "opus", "aac", "flac", "wav", "pcm"] | None = None
    speed: Annotated[float | None, Field(ge=0.25, le=4.0)] = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
