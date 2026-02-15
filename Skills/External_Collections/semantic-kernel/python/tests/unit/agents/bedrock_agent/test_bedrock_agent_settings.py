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

import pytest
from pydantic import ValidationError

from semantic_kernel.agents.bedrock.bedrock_agent_settings import BedrockAgentSettings


def test_bedrock_agent_settings_from_env_vars(bedrock_agent_unit_test_env):
    """Test loading BedrockAgentSettings from environment variables."""
    settings = BedrockAgentSettings(env_file_path="fake_path")

    assert settings.agent_resource_role_arn == bedrock_agent_unit_test_env["BEDROCK_AGENT_AGENT_RESOURCE_ROLE_ARN"]
    assert settings.foundation_model == bedrock_agent_unit_test_env["BEDROCK_AGENT_FOUNDATION_MODEL"]


@pytest.mark.parametrize(
    "exclude_list",
    [
        ["BEDROCK_AGENT_AGENT_RESOURCE_ROLE_ARN"],
        ["BEDROCK_AGENT_FOUNDATION_MODEL"],
    ],
    indirect=True,
)
def test_bedrock_agent_settings_from_env_vars_missing_required(bedrock_agent_unit_test_env):
    """Test loading BedrockAgentSettings from environment variables with missing required fields."""
    with pytest.raises(ValidationError):
        BedrockAgentSettings(env_file_path="fake_path")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
