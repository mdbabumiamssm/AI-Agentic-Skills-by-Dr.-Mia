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

import sys
from typing import Any

from openai import _legacy_response

if sys.version_info >= (3, 12):
    from typing import override  # pragma: no cover
else:
    from typing_extensions import override  # pragma: no cover

from semantic_kernel.connectors.ai.open_ai.prompt_execution_settings.open_ai_text_to_audio_execution_settings import (
    OpenAITextToAudioExecutionSettings,
)
from semantic_kernel.connectors.ai.open_ai.services.open_ai_handler import OpenAIHandler
from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings
from semantic_kernel.connectors.ai.text_to_audio_client_base import TextToAudioClientBase
from semantic_kernel.contents.audio_content import AudioContent


class OpenAITextToAudioBase(OpenAIHandler, TextToAudioClientBase):
    """OpenAI text to audio client base class."""

    @override
    async def get_audio_contents(
        self,
        text: str,
        settings: PromptExecutionSettings | None = None,
        **kwargs: Any,
    ) -> list[AudioContent]:
        if not settings:
            settings = OpenAITextToAudioExecutionSettings(ai_model_id=self.ai_model_id)
        else:
            if not isinstance(settings, OpenAITextToAudioExecutionSettings):
                settings = self.get_prompt_execution_settings_from_settings(settings)

        assert isinstance(settings, OpenAITextToAudioExecutionSettings)  # nosec

        if settings.ai_model_id is None:
            settings.ai_model_id = self.ai_model_id
        settings.input = text

        response = await self._send_request(settings)
        assert isinstance(response, _legacy_response.HttpxBinaryResponseContent)  # nosec

        return [
            AudioContent(
                ai_model_id=settings.ai_model_id,
                data=response.read(),
                data_format="base64",
            )
        ]

    def get_prompt_execution_settings_class(self) -> type[PromptExecutionSettings]:
        """Get the request settings class."""
        return OpenAITextToAudioExecutionSettings

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
