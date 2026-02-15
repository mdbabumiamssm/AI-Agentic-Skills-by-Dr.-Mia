# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

from importlib import metadata
from importlib.metadata import PackageMetadata

from squidpy import datasets, experimental, gr, im, pl, read, tl

try:
    md: PackageMetadata = metadata.metadata(__name__)
    __version__ = md["Version"] if "Version" in md else ""
    __author__ = md["Author"] if "Author" in md else ""
    __maintainer__ = md["Maintainer-email"] if "Maintainer-email" in md else ""
except ImportError:
    md = None  # type: ignore[assignment]

del metadata, md

__all__ = ["datasets", "experimental", "gr", "im", "pl", "read", "tl"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
