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


from unittest.mock import patch

import pytest

from semantic_kernel.agents.chat_completion.chat_completion_agent import ChatCompletionAgent, ChatHistoryAgentThread
from semantic_kernel.exceptions.kernel_exceptions import KernelServiceNotFoundError


@patch("semantic_kernel.utils.telemetry.agent_diagnostics.decorators.tracer")
async def test_chat_completion_agent_get_response(
    mock_tracer,
    chat_history,
    model_diagnostics_unit_test_env,
):
    # Arrange
    chat_completion_agent = ChatCompletionAgent()
    thread = ChatHistoryAgentThread(chat_history=chat_history)
    # Act
    with pytest.raises(KernelServiceNotFoundError):
        await chat_completion_agent.get_response(messages="test", thread=thread)
    # Assert
    mock_tracer.start_as_current_span.assert_called_once()
    args, _ = mock_tracer.start_as_current_span.call_args
    assert args[0] == f"invoke_agent {chat_completion_agent.name}"


@patch("semantic_kernel.utils.telemetry.agent_diagnostics.decorators.tracer")
async def test_chat_completion_agent_invoke(
    mock_tracer,
    chat_history,
    model_diagnostics_unit_test_env,
):
    # Arrange
    chat_completion_agent = ChatCompletionAgent()
    thread = ChatHistoryAgentThread(chat_history=chat_history)
    # Act
    with pytest.raises(KernelServiceNotFoundError):
        async for _ in chat_completion_agent.invoke(messages="test", thread=thread):
            pass
    # Assert
    mock_tracer.start_as_current_span.assert_called_once()
    args, _ = mock_tracer.start_as_current_span.call_args
    assert args[0] == f"invoke_agent {chat_completion_agent.name}"


@patch("semantic_kernel.utils.telemetry.agent_diagnostics.decorators.tracer")
async def test_chat_completion_agent_invoke_stream(
    mock_tracer,
    chat_history,
    model_diagnostics_unit_test_env,
):
    # Arrange
    chat_completion_agent = ChatCompletionAgent()
    thread = ChatHistoryAgentThread(chat_history=chat_history)
    # Act
    with pytest.raises(KernelServiceNotFoundError):
        async for _ in chat_completion_agent.invoke_stream(messages="test", thread=thread):
            pass
    # Assert
    mock_tracer.start_as_current_span.assert_called_once()
    args, _ = mock_tracer.start_as_current_span.call_args
    assert args[0] == f"invoke_agent {chat_completion_agent.name}"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
