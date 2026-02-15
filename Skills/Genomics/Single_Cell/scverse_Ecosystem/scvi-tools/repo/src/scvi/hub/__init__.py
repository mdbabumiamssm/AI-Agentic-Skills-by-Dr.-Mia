# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from scvi.utils import error_on_missing_dependencies

error_on_missing_dependencies("huggingface_hub")


from ._metadata import HubMetadata, HubModelCardHelper  # noqa: E402
from ._model import HubModel  # noqa: E402

__all__ = [
    "HubModel",
    "HubMetadata",
    "HubModelCardHelper",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
