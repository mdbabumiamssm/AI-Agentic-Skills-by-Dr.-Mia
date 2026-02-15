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
from typing import Any

from pydantic import Field

from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings

logger = logging.getLogger(__name__)


class OpenAIAudioToTextExecutionSettings(PromptExecutionSettings):
    """Request settings for OpenAI audio to text services."""

    ai_model_id: str | None = Field(default=None, serialization_alias="model")
    filename: str | None = Field(
        default=None,
        description="Do not set this manually. It is set by the service based on the audio content.",
    )
    language: str | None = None
    prompt: str | None = None
    response_format: str | None = None
    temperature: float | None = None

    def prepare_settings_dict(self, **kwargs) -> dict[str, Any]:
        """Prepare the settings dictionary for the OpenAI API."""
        settings_dict = super().prepare_settings_dict(**kwargs)

        # Remove the file name since it will be open as a file object
        settings_dict.pop("filename", None)

        return settings_dict

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
