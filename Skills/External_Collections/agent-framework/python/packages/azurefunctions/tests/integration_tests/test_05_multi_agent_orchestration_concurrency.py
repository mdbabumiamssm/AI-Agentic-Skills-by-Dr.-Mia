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
"""
Integration Tests for MultiAgent Concurrency Sample

Tests the multi-agent concurrency sample for parallel agent execution.

The function app is automatically started by the test fixture.

Prerequisites:
- Azure OpenAI credentials configured (see packages/azurefunctions/tests/integration_tests/.env.example)
- Azurite running for durable orchestrations (or Azure Storage account configured)

Usage:
    # Start Azurite (if not already running)
    azurite &

    # Run tests
    uv run pytest packages/azurefunctions/tests/integration_tests/test_05_multi_agent_orchestration_concurrency.py -v
"""

import pytest

from .testutils import SampleTestHelper, skip_if_azure_functions_integration_tests_disabled

# Module-level markers - applied to all tests in this file
pytestmark = [
    pytest.mark.orchestration,
    pytest.mark.sample("05_multi_agent_orchestration_concurrency"),
    pytest.mark.usefixtures("function_app_for_test"),
    skip_if_azure_functions_integration_tests_disabled,
]


class TestSampleMultiAgentConcurrency:
    """Tests for 05_multi_agent_orchestration_concurrency sample."""

    def test_concurrent_agents(self, base_url: str) -> None:
        """Test multiple agents running concurrently."""
        # Start orchestration
        response = SampleTestHelper.post_text(f"{base_url}/api/multiagent/run", "What is temperature?")
        assert response.status_code == 202
        data = response.json()
        assert "instanceId" in data
        assert "statusQueryGetUri" in data

        # Wait for completion
        status = SampleTestHelper.wait_for_orchestration(data["statusQueryGetUri"])
        assert status["runtimeStatus"] == "Completed"
        output = status["output"]
        assert "physicist" in output
        assert "chemist" in output


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
