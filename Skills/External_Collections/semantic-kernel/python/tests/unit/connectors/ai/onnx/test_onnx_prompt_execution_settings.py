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

from semantic_kernel.connectors.ai.onnx.onnx_gen_ai_prompt_execution_settings import (
    OnnxGenAIPromptExecutionSettings,
)
from semantic_kernel.connectors.ai.prompt_execution_settings import PromptExecutionSettings


def test_default_onnx_chat_prompt_execution_settings():
    settings = OnnxGenAIPromptExecutionSettings()
    assert settings.temperature is None
    assert settings.top_p is None


def test_custom_onnx_chat_prompt_execution_settings():
    settings = OnnxGenAIPromptExecutionSettings(
        temperature=0.5,
        top_p=0.5,
        max_length=128,
    )
    assert settings.temperature == 0.5
    assert settings.top_p == 0.5
    assert settings.max_length == 128


def test_onnx_chat_prompt_execution_settings_from_default_completion_config():
    settings = PromptExecutionSettings(service_id="test_service")
    chat_settings = OnnxGenAIPromptExecutionSettings.from_prompt_execution_settings(settings)
    assert chat_settings.service_id == "test_service"
    assert chat_settings.temperature is None
    assert chat_settings.top_p is None


def test_onnx_chat_prompt_execution_settings_from_onnx_prompt_execution_settings():
    chat_settings = OnnxGenAIPromptExecutionSettings(service_id="test_service", temperature=1.0)
    new_settings = OnnxGenAIPromptExecutionSettings(service_id="test_2", temperature=0.0)
    chat_settings.update_from_prompt_execution_settings(new_settings)
    assert chat_settings.service_id == "test_2"
    assert chat_settings.temperature == 0.0


def test_onnx_chat_prompt_execution_settings_from_custom_completion_config():
    settings = PromptExecutionSettings(
        service_id="test_service",
        extension_data={
            "temperature": 0.5,
            "top_p": 0.5,
            "max_length": 128,
        },
    )
    chat_settings = OnnxGenAIPromptExecutionSettings.from_prompt_execution_settings(settings)
    assert chat_settings.temperature == 0.5
    assert chat_settings.top_p == 0.5
    assert chat_settings.max_length == 128


def test_create_options():
    settings = OnnxGenAIPromptExecutionSettings(
        service_id="test_service",
        extension_data={
            "temperature": 0.5,
            "top_p": 0.5,
            "max_length": 128,
        },
    )
    options = settings.prepare_settings_dict()
    assert options["temperature"] == 0.5
    assert options["top_p"] == 0.5
    assert options["max_length"] == 128


def test_create_options_with_wrong_parameter():
    with pytest.raises(ValidationError):
        OnnxGenAIPromptExecutionSettings(
            service_id="test_service",
            function_choice_behavior="auto",
            extension_data={
                "temperature": 10.0,
                "top_p": 0.5,
                "max_length": 128,
            },
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
