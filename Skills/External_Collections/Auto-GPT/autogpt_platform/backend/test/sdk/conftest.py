# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""
Configuration for SDK tests.

This conftest.py file provides basic test setup for SDK unit tests
without requiring the full server infrastructure.
"""

from unittest.mock import MagicMock

import pytest


@pytest.fixture(scope="session")
def server():
    """Mock server fixture for SDK tests."""
    mock_server = MagicMock()
    mock_server.agent_server = MagicMock()
    mock_server.agent_server.test_create_graph = MagicMock()
    return mock_server


@pytest.fixture(autouse=True)
def reset_registry():
    """Reset the AutoRegistry before each test."""
    from backend.sdk.registry import AutoRegistry

    AutoRegistry.clear()
    yield
    AutoRegistry.clear()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
