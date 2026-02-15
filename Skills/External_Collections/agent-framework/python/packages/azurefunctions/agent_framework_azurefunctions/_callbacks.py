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

"""Callback interfaces for Durable Agent executions.

This module enables callers of AgentFunctionApp to supply streaming and final-response callbacks that are
invoked during durable entity execution.
"""

from dataclasses import dataclass
from typing import Protocol

from agent_framework import AgentRunResponse, AgentRunResponseUpdate


@dataclass(frozen=True)
class AgentCallbackContext:
    """Context supplied to callback invocations."""

    agent_name: str
    correlation_id: str
    thread_id: str | None = None
    request_message: str | None = None


class AgentResponseCallbackProtocol(Protocol):
    """Protocol describing the callbacks invoked during agent execution."""

    async def on_streaming_response_update(
        self,
        update: AgentRunResponseUpdate,
        context: AgentCallbackContext,
    ) -> None:
        """Handle a streaming response update emitted by the agent."""

    async def on_agent_response(
        self,
        response: AgentRunResponse,
        context: AgentCallbackContext,
    ) -> None:
        """Handle the final agent response."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
