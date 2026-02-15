# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from backend.server.v2.library.model import LibraryAgentPreset

from .graph import NodeModel
from .integrations import Webhook  # noqa: F401

# Resolve Webhook forward references
NodeModel.model_rebuild()
LibraryAgentPreset.model_rebuild()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
