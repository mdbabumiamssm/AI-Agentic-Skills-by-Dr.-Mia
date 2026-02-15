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
"""Stream responses from Copilot Studio agents in SK and AF."""

import asyncio


async def run_semantic_kernel() -> None:
    from semantic_kernel.agents import CopilotStudioAgent

    agent = CopilotStudioAgent(
        name="TourGuide",
        instructions="Provide travel recommendations in short bursts.",
    )
    # SK streaming yields chunks with message metadata.
    print("[SK][stream]", end=" ")
    async for chunk in agent.invoke_stream("Plan a day in Copenhagen for foodies."):
        if chunk.message:
            print(chunk.message.content, end="", flush=True)
    print()


async def run_agent_framework() -> None:
    from agent_framework.microsoft import CopilotStudioAgent

    agent = CopilotStudioAgent(
        name="TourGuide",
        instructions="Provide travel recommendations in short bursts.",
    )
    # AF streaming provides incremental AgentRunResponseUpdate objects.
    print("[AF][stream]", end=" ")
    async for update in agent.run_stream("Plan a day in Copenhagen for foodies."):
        if update.text:
            print(update.text, end="", flush=True)
    print()


async def main() -> None:
    await run_semantic_kernel()
    await run_agent_framework()


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
