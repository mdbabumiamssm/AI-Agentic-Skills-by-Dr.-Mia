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
from typing import ClassVar

from microsoft_agents.copilotstudio.client import (
    AgentType,
    PowerPlatformCloud,
)
from pydantic import Field, SecretStr

from semantic_kernel.kernel_pydantic import KernelBaseSettings
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class CopilotStudioAgentAuthMode(str, Enum):
    """The Copilot Studio agent authentication mode."""

    INTERACTIVE = "interactive"  # user authentication
    SERVICE = "service"  # client-credentials (app secret/cert)


@experimental
class CopilotStudioAgentSettings(KernelBaseSettings):
    """Copilot Studio Agent settings currently used by the CopilotStudioAgent."""

    env_prefix: ClassVar[str] = "COPILOT_STUDIO_AGENT_"

    app_client_id: str | None = None
    tenant_id: str | None = None
    environment_id: str | None = None
    agent_identifier: str | None = None
    cloud: PowerPlatformCloud = Field(default=PowerPlatformCloud.UNKNOWN)
    copilot_agent_type: AgentType = Field(default=AgentType.PUBLISHED)
    custom_power_platform_cloud: str | None = None
    client_secret: SecretStr | None = None
    client_certificate: str | None = None
    user_assertion: str | None = None
    auth_mode: CopilotStudioAgentAuthMode = Field(default=CopilotStudioAgentAuthMode.INTERACTIVE)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
