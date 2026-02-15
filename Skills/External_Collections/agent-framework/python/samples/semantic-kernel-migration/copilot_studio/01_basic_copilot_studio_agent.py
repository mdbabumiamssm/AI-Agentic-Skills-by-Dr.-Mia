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
"""Call a Copilot Studio agent with SK and Agent Framework."""

import asyncio


async def run_semantic_kernel() -> None:
    from semantic_kernel.agents import CopilotStudioAgent

    # SK agent talks to the configured Copilot Studio bot directly.
    agent = CopilotStudioAgent(
        name="PhysicsAgent",
        instructions="Answer physics questions concisely.",
    )
    response = await agent.get_response("Why is the sky blue?")
    print("[SK]", response.message.content)


async def run_agent_framework() -> None:
    from agent_framework.microsoft import CopilotStudioAgent

    # AF exposes an equivalent CopilotStudioAgent wrapper.
    agent = CopilotStudioAgent(
        name="PhysicsAgent",
        instructions="Answer physics questions concisely.",
    )
    reply = await agent.run("Why is the sky blue?")
    print("[AF]", reply.text)


async def main() -> None:
    await run_semantic_kernel()
    await run_agent_framework()


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
