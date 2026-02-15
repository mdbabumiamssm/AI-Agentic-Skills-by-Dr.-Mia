# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

""" GTTS Voice. """
from __future__ import annotations

import os

import gtts
from playsound import playsound

from .base import VoiceBase


class GTTSVoice(VoiceBase):
    """GTTS Voice."""

    def _setup(self) -> None:
        pass

    def _speech(self, text: str, voice_id: int = 0) -> bool:
        """Play the given text."""
        tts = gtts.gTTS(text)
        tts.save("speech.mp3")
        playsound("speech.mp3", True)
        os.remove("speech.mp3")
        return True

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
