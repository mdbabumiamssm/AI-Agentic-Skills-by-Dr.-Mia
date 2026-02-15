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


class FilterTypes(str, Enum):
    """Enum for the filter types."""

    FUNCTION_INVOCATION = "function_invocation"
    PROMPT_RENDERING = "prompt_rendering"
    AUTO_FUNCTION_INVOCATION = "auto_function_invocation"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
