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

from enum import Enum


class OpenAIModelTypes(Enum):
    """OpenAI model types, can be text, chat or embedding."""

    TEXT = "text"
    CHAT = "chat"
    EMBEDDING = "embedding"
    TEXT_TO_IMAGE = "text-to-image"
    AUDIO_TO_TEXT = "audio-to-text"
    TEXT_TO_AUDIO = "text-to-audio"
    REALTIME = "realtime"
    RESPONSE = "response"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
