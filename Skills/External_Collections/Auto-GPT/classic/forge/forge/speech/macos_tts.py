# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

""" MacOS TTS Voice. """
from __future__ import annotations

import subprocess

from .base import VoiceBase


class MacOSTTS(VoiceBase):
    """MacOS TTS Voice."""

    def _setup(self) -> None:
        pass

    def _speech(self, text: str, voice_id: int = 0) -> bool:
        """Play the given text."""
        if voice_id == 0:
            subprocess.run(["say", text], shell=False)
        elif voice_id == 1:
            subprocess.run(["say", "-v", "Ava (Premium)", text], shell=False)
        else:
            subprocess.run(["say", "-v", "Samantha", text], shell=False)
        return True

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
