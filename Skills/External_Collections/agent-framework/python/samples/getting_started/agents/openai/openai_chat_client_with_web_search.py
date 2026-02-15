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

import asyncio

from agent_framework import ChatAgent, HostedWebSearchTool
from agent_framework.openai import OpenAIChatClient

"""
OpenAI Chat Client with Web Search Example

This sample demonstrates using HostedWebSearchTool with OpenAI Chat Client
for real-time information retrieval and current data access.
"""


async def main() -> None:
    # Test that the agent will use the web search tool with location
    additional_properties = {
        "user_location": {
            "country": "US",
            "city": "Seattle",
        }
    }

    agent = ChatAgent(
        chat_client=OpenAIChatClient(model_id="gpt-4o-search-preview"),
        instructions="You are a helpful assistant that can search the web for current information.",
        tools=[HostedWebSearchTool(additional_properties=additional_properties)],
    )

    message = "What is the current weather? Do not ask for my current location."
    stream = False
    print(f"User: {message}")

    if stream:
        print("Assistant: ", end="")
        async for chunk in agent.run_stream(message):
            if chunk.text:
                print(chunk.text, end="")
        print("")
    else:
        response = await agent.run(message)
        print(f"Assistant: {response}")


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
