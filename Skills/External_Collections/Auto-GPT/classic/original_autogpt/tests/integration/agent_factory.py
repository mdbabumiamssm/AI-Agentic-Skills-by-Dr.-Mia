# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pathlib import Path

import pytest
from forge.config.ai_profile import AIProfile
from forge.file_storage import FileStorageBackendName, get_storage
from forge.llm.providers import MultiProvider

from autogpt.agents.agent import Agent, AgentConfiguration, AgentSettings
from autogpt.app.config import AppConfig


@pytest.fixture
def dummy_agent(config: AppConfig, llm_provider: MultiProvider):
    ai_profile = AIProfile(
        ai_name="Dummy Agent",
        ai_role="Dummy Role",
        ai_goals=[
            "Dummy Task",
        ],
    )

    agent_settings = AgentSettings(
        name=Agent.default_settings.name,
        description=Agent.default_settings.description,
        ai_profile=ai_profile,
        config=AgentConfiguration(
            fast_llm=config.fast_llm,
            smart_llm=config.smart_llm,
            use_functions_api=config.openai_functions,
        ),
        history=Agent.default_settings.history.model_copy(deep=True),
    )

    local = config.file_storage_backend == FileStorageBackendName.LOCAL
    restrict_to_root = not local or config.restrict_to_workspace
    file_storage = get_storage(
        config.file_storage_backend,
        root_path=Path("data"),
        restrict_to_root=restrict_to_root,
    )
    file_storage.initialize()

    agent = Agent(
        settings=agent_settings,
        llm_provider=llm_provider,
        file_storage=file_storage,
        app_config=config,
    )

    return agent

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
