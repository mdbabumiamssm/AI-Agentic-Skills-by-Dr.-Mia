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

from typing import Annotated, Any

import anyio
from agent_framework.openai import OpenAIResponsesClient

"""
This sample demonstrates how to expose an Agent as an MCP server.

To run this sample, set up your MCP host (like Claude Desktop or VSCode Github Copilot Agents)
with the following configuration:
```json
{
    "servers": {
        "agent-framework": {
            "command": "uv",
            "args": [
                "--directory=<path to project>/agent-framework/python/samples/getting_started/mcp",
                "run",
                "agent_as_mcp_server.py"
            ],
            "env": {
                "OPENAI_API_KEY": "<OpenAI API key>",
                "OPENAI_RESPONSES_MODEL_ID": "<OpenAI Responses model ID>",
            }
        }
    }
}
```
"""


def get_specials() -> Annotated[str, "Returns the specials from the menu."]:
    return """
        Special Soup: Clam Chowder
        Special Salad: Cobb Salad
        Special Drink: Chai Tea
        """


def get_item_price(
    menu_item: Annotated[str, "The name of the menu item."],
) -> Annotated[str, "Returns the price of the menu item."]:
    return "$9.99"


async def run() -> None:
    # Define an agent
    # Agent's name and description provide better context for AI model
    agent = OpenAIResponsesClient().create_agent(
        name="RestaurantAgent",
        description="Answer questions about the menu.",
        tools=[get_specials, get_item_price],
    )

    # Expose the agent as an MCP server
    server = agent.as_mcp_server()

    # Run server
    from mcp.server.stdio import stdio_server

    async def handle_stdin(stdin: Any | None = None, stdout: Any | None = None) -> None:
        async with stdio_server() as (read_stream, write_stream):
            await server.run(read_stream, write_stream, server.create_initialization_options())

    await handle_stdin()


if __name__ == "__main__":
    anyio.run(run)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
