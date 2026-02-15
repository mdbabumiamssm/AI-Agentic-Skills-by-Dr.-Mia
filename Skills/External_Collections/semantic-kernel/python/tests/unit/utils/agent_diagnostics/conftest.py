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

import semantic_kernel
from semantic_kernel.utils.telemetry.model_diagnostics.model_diagnostics_settings import ModelDiagnosticSettings


@pytest.fixture()
def model_diagnostics_unit_test_env(monkeypatch):
    """Fixture to set environment variables for Model Diagnostics Unit Tests."""
    env_vars = {
        "SEMANTICKERNEL_EXPERIMENTAL_GENAI_ENABLE_OTEL_DIAGNOSTICS": "true",
        "SEMANTICKERNEL_EXPERIMENTAL_GENAI_ENABLE_OTEL_DIAGNOSTICS_SENSITIVE": "true",
    }

    for key, value in env_vars.items():
        monkeypatch.setenv(key, value)

    # Need to reload the settings to pick up the new environment variables since the
    # settings are loaded at import time and this fixture is called after the import
    semantic_kernel.utils.telemetry.agent_diagnostics.decorators.MODEL_DIAGNOSTICS_SETTINGS = ModelDiagnosticSettings()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
