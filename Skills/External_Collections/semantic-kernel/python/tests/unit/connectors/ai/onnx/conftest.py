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
from pytest import fixture


@fixture()
def onnx_unit_test_env(monkeypatch, exclude_list, override_env_param_dict):
    """Fixture to set environment variables for OnnxGenAISettings."""
    if exclude_list is None:
        exclude_list = []

    if override_env_param_dict is None:
        override_env_param_dict = {}

    env_vars = {
        "ONNX_GEN_AI_CHAT_MODEL_FOLDER": "test",
        "ONNX_GEN_AI_TEXT_MODEL_FOLDER": "test",
    }

    env_vars.update(override_env_param_dict)

    for key, value in env_vars.items():
        if key not in exclude_list:
            monkeypatch.setenv(key, value)
        else:
            monkeypatch.delenv(key, raising=False)

    return env_vars


gen_ai_config = {"model": {"test": "test"}}

gen_ai_config_vision = {"model": {"vision": "test"}}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
