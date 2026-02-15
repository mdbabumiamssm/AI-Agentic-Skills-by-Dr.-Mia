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

import pytest
from azure.identity import AzureCliCredential

from semantic_kernel.connectors.ai.open_ai.services.azure_text_to_image import AzureTextToImage
from semantic_kernel.connectors.ai.open_ai.services.open_ai_text_to_image import OpenAITextToImage
from semantic_kernel.connectors.ai.text_to_image_client_base import TextToImageClientBase


class TextToImageTestBase:
    """Base class for testing text-to-image services."""

    @pytest.fixture(scope="module")
    def services(self) -> dict[str, TextToImageClientBase]:
        """Return text-to-image services."""
        return {
            "openai": OpenAITextToImage(),
            "azure_openai": AzureTextToImage(credential=AzureCliCredential()),
        }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
