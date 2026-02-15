# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from forge.file_storage import FileStorageBackendName
from forge.models.config import SystemSettings, UserConfigurable
from forge.speech.say import TTSConfig


class BaseConfig(SystemSettings):
    name: str = "Base configuration"
    description: str = "Default configuration for forge agent."

    # TTS configuration
    tts_config: TTSConfig = TTSConfig()

    # File storage
    file_storage_backend: FileStorageBackendName = UserConfigurable(
        default=FileStorageBackendName.LOCAL, from_env="FILE_STORAGE_BACKEND"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
