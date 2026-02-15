# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import enum
from pathlib import Path

from .base import FileStorage


class FileStorageBackendName(str, enum.Enum):
    LOCAL = "local"
    GCS = "gcs"
    S3 = "s3"


def get_storage(
    backend: FileStorageBackendName,
    root_path: Path = Path("."),
    restrict_to_root: bool = True,
) -> FileStorage:
    match backend:
        case FileStorageBackendName.LOCAL:
            from .local import FileStorageConfiguration, LocalFileStorage

            config = FileStorageConfiguration.from_env()
            config.root = root_path
            config.restrict_to_root = restrict_to_root
            return LocalFileStorage(config)
        case FileStorageBackendName.S3:
            from .s3 import S3FileStorage, S3FileStorageConfiguration

            config = S3FileStorageConfiguration.from_env()
            config.root = root_path
            return S3FileStorage(config)
        case FileStorageBackendName.GCS:
            from .gcs import GCSFileStorage, GCSFileStorageConfiguration

            config = GCSFileStorageConfiguration.from_env()
            config.root = root_path
            return GCSFileStorage(config)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
