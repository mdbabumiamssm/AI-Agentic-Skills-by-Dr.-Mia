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

"""Simple agentic chat example (Feature 1: Agentic Chat)."""

from agent_framework import ChatAgent, ChatClientProtocol


def simple_agent(chat_client: ChatClientProtocol) -> ChatAgent:
    """Create a simple chat agent.

    Args:
        chat_client: The chat client to use for the agent

    Returns:
        A configured ChatAgent instance
    """
    return ChatAgent(
        name="simple_chat_agent",
        instructions="You are a helpful assistant. Be concise and friendly.",
        chat_client=chat_client,
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
