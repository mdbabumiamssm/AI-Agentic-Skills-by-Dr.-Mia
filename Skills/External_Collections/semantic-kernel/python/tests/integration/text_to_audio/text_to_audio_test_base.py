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

import os

import pytest
from azure.identity import AzureCliCredential

from semantic_kernel.connectors.ai.open_ai import AzureTextToAudio, OpenAITextToAudio
from semantic_kernel.connectors.ai.text_to_audio_client_base import TextToAudioClientBase
from tests.utils import is_service_setup_for_testing

# TTS model on Azure model is not available in regions at which we have chat completion models.
# Therefore, we need to use a different endpoint for testing.
azure_setup = is_service_setup_for_testing(["AZURE_OPENAI_TEXT_TO_AUDIO_ENDPOINT"])


class TextToAudioTestBase:
    """Base class for testing text-to-audio services."""

    @pytest.fixture(scope="module")
    def services(self) -> dict[str, TextToAudioClientBase]:
        """Return text-to-audio services."""
        return {
            "openai": OpenAITextToAudio(),
            "azure_openai": AzureTextToAudio(
                endpoint=os.environ["AZURE_OPENAI_TEXT_TO_AUDIO_ENDPOINT"], credential=AzureCliCredential()
            )
            if azure_setup
            else None,
        }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
