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

from agent_framework import HostedMCPTool
from agent_framework.azure import AzureOpenAIChatClient
from azure.ai.agentserver.agentframework import from_agent_framework  # pyright: ignore[reportUnknownVariableType]
from azure.identity import DefaultAzureCredential


def main():
    # Create an Agent using the Azure OpenAI Chat Client with a MCP Tool that connects to Microsoft Learn MCP
    agent = AzureOpenAIChatClient(credential=DefaultAzureCredential()).create_agent(
        name="DocsAgent",
        instructions="You are a helpful assistant that can help with microsoft documentation questions.",
        tools=HostedMCPTool(
            name="Microsoft Learn MCP",
            url="https://learn.microsoft.com/api/mcp",
        ),
    )

    # Run the agent as a hosted agent
    from_agent_framework(agent).run()


if __name__ == "__main__":
    main()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
