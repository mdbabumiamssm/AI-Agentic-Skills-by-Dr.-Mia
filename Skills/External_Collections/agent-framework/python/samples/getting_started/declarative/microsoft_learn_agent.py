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
from pathlib import Path

from agent_framework.declarative import AgentFactory
from azure.identity.aio import AzureCliCredential


async def main():
    """Create an agent from a declarative yaml specification and run it."""
    # get the path
    current_path = Path(__file__).parent
    yaml_path = current_path.parent.parent.parent.parent / "agent-samples" / "foundry" / "MicrosoftLearnAgent.yaml"

    # create the agent from the yaml
    async with (
        AzureCliCredential() as credential,
        AgentFactory(client_kwargs={"credential": credential}).create_agent_from_yaml_path(yaml_path) as agent,
    ):
        response = await agent.run("How do I create a storage account with private endpoint using bicep?")
        print("Agent response:", response.text)


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
