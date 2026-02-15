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
import os

from agent_framework.azure import AzureAIClient
from azure.identity.aio import AzureCliCredential

"""
Azure AI Agent with Microsoft Fabric Example

This sample demonstrates usage of AzureAIClient with Microsoft Fabric
to query Fabric data sources and provide responses based on data analysis.

Prerequisites:
1. Set AZURE_AI_PROJECT_ENDPOINT and AZURE_AI_MODEL_DEPLOYMENT_NAME environment variables.
2. Ensure you have a Microsoft Fabric connection configured in your Azure AI project
   and set FABRIC_PROJECT_CONNECTION_ID environment variable.
"""


async def main() -> None:
    async with (
        AzureCliCredential() as credential,
        AzureAIClient(credential=credential).create_agent(
            name="MyFabricAgent",
            instructions="You are a helpful assistant.",
            tools={
                "type": "fabric_dataagent_preview",
                "fabric_dataagent_preview": {
                    "project_connections": [
                        {
                            "project_connection_id": os.environ["FABRIC_PROJECT_CONNECTION_ID"],
                        }
                    ]
                },
            },
        ) as agent,
    ):
        query = "Tell me about sales records"
        print(f"User: {query}")
        result = await agent.run(query)
        print(f"Result: {result}\n")


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
