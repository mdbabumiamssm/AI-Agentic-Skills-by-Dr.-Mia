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

from dataclasses import dataclass

from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
@dataclass
class RestApiOAuthFlow:
    """Represents the OAuth flow used by the REST API."""

    authorization_url: str
    token_url: str
    scopes: dict[str, str]
    refresh_url: str | None = None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
