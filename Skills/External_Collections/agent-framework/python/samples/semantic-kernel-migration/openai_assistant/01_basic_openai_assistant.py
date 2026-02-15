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
"""Create an OpenAI Assistant using SK and Agent Framework."""

import asyncio
import os

ASSISTANT_MODEL = os.environ.get("OPENAI_ASSISTANT_MODEL", "gpt-4o-mini")


async def run_semantic_kernel() -> None:
    from semantic_kernel.agents import AssistantAgentThread, OpenAIAssistantAgent

    client = OpenAIAssistantAgent.create_client()
    # Provision the assistant on the OpenAI Assistants service.
    definition = await client.beta.assistants.create(
        model=ASSISTANT_MODEL,
        name="Helper",
        instructions="Answer questions in one concise paragraph.",
    )
    agent = OpenAIAssistantAgent(client=client, definition=definition)

    thread: AssistantAgentThread | None = None
    response = await agent.get_response("What is the capital of Denmark?", thread=thread)
    thread = response.thread
    print("[SK]", response.message.content)
    if thread is not None:
        print("[SK][thread-id]", thread.id)


async def run_agent_framework() -> None:
    from agent_framework.openai import OpenAIAssistantsClient

    assistants_client = OpenAIAssistantsClient()
    # AF wraps the assistant lifecycle with an async context manager.
    async with assistants_client.create_agent(
        name="Helper",
        instructions="Answer questions in one concise paragraph.",
        model=ASSISTANT_MODEL,
    ) as assistant_agent:
        reply = await assistant_agent.run("What is the capital of Denmark?")
        print("[AF]", reply.text)
        follow_up = await assistant_agent.run(
            "How many residents live there?",
            thread=assistant_agent.get_new_thread(),
        )
        print("[AF][follow-up]", follow_up.text)


async def main() -> None:
    await run_semantic_kernel()
    await run_agent_framework()


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
