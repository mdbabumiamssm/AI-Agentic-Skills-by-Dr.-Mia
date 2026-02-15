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

import re

from semantic_kernel.utils.feature_stage_decorator import experimental

_AGENT_TYPE_REGEX = re.compile(r"^[\w\-\.]+\Z")


@experimental
def is_valid_agent_type(value: str) -> bool:
    """Check if the agent type is valid."""
    return bool(_AGENT_TYPE_REGEX.match(value))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
