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

from typing import ClassVar

from pydantic import SecretStr

from semantic_kernel.kernel_pydantic import KernelBaseSettings


class CrewAISettings(KernelBaseSettings):
    """The Crew.AI settings.

    Required:
    - endpoint: str - The API endpoint.
    """

    env_prefix: ClassVar[str] = "CREW_AI_"

    endpoint: str
    auth_token: SecretStr
    polling_interval: float = 1.0
    polling_timeout: float = 30.0

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
