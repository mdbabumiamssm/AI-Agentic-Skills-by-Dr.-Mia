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

from unittest.mock import AsyncMock, MagicMock

import pytest
from azure.ai.agents.models import Agent as AzureAIAgentModel
from azure.ai.projects.aio import AIProjectClient


@pytest.fixture
def ai_project_client() -> AsyncMock:
    client = AsyncMock(spec=AIProjectClient)

    agents_mock = MagicMock()
    client.agents = agents_mock

    return client


@pytest.fixture
def ai_agent_definition() -> AsyncMock:
    definition = AsyncMock(spec=AzureAIAgentModel)
    definition.id = "agent123"
    definition.name = "agentName"
    definition.description = "desc"
    definition.instructions = "test agent"

    return definition

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
