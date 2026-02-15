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
Integration Tests for Multi-Agent Sample

Tests the multi-agent sample with different agent endpoints.

The function app is automatically started by the test fixture.

Prerequisites:
- Azure OpenAI credentials configured (see packages/azurefunctions/tests/integration_tests/.env.example)
- Azurite or Azure Storage account configured

Usage:
    uv run pytest packages/azurefunctions/tests/integration_tests/test_02_multi_agent.py -v
"""

import pytest

from .testutils import SampleTestHelper, skip_if_azure_functions_integration_tests_disabled

# Module-level markers - applied to all tests in this file
pytestmark = [
    pytest.mark.sample("02_multi_agent"),
    pytest.mark.usefixtures("function_app_for_test"),
    skip_if_azure_functions_integration_tests_disabled,
]


class TestSampleMultiAgent:
    """Tests for 02_multi_agent sample."""

    @pytest.fixture(autouse=True)
    def _set_agent_urls(self, base_url: str) -> None:
        """Configure base URLs for Weather and Math agents."""
        self.weather_base_url = f"{base_url}/api/agents/WeatherAgent"
        self.math_base_url = f"{base_url}/api/agents/MathAgent"

    def test_weather_agent(self) -> None:
        """Test WeatherAgent endpoint."""
        response = SampleTestHelper.post_json(
            f"{self.weather_base_url}/run",
            {"message": "What is the weather in Seattle?"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "response" in data

    def test_math_agent(self) -> None:
        """Test MathAgent endpoint."""
        response = SampleTestHelper.post_json(
            f"{self.math_base_url}/run",
            {"message": "Calculate a 20% tip on a $50 bill", "wait_for_response": False},
        )
        assert response.status_code == 202
        data = response.json()

        assert data["status"] == "accepted"
        assert "correlation_id" in data
        assert "thread_id" in data


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
