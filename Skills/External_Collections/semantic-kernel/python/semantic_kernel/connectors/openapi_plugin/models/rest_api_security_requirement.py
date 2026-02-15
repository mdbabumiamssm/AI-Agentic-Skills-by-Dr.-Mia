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

from semantic_kernel.connectors.openapi_plugin.models.rest_api_security_scheme import RestApiSecurityScheme
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class RestApiSecurityRequirement(dict[RestApiSecurityScheme, list[str]]):
    """Represents the security requirements used by the REST API."""

    def __init__(self, dictionary: dict[RestApiSecurityScheme, list[str]]):
        """Initializes a new instance of the RestApiSecurityRequirement class."""
        super().__init__(dictionary)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
