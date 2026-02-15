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

from agent_framework.azure import AzureAIClient
from azure.identity.aio import AzureCliCredential
from pydantic import BaseModel, ConfigDict

"""
Azure AI Agent Response Format Example

This sample demonstrates basic usage of AzureAIClient with response format,
also known as structured outputs.
"""


class ReleaseBrief(BaseModel):
    feature: str
    benefit: str
    launch_date: str
    model_config = ConfigDict(extra="forbid")


async def main() -> None:
    """Example of using response_format property."""

    # Since no Agent ID is provided, the agent will be automatically created.
    # For authentication, run `az login` command in terminal or replace AzureCliCredential with preferred
    # authentication option.
    async with (
        AzureCliCredential() as credential,
        AzureAIClient(credential=credential).create_agent(
            name="ProductMarketerAgent",
            instructions="Return launch briefs as structured JSON.",
        ) as agent,
    ):
        query = "Draft a launch brief for the Contoso Note app."
        print(f"User: {query}")
        result = await agent.run(
            query,
            # Specify type to use as response
            response_format=ReleaseBrief,
        )

        if isinstance(result.value, ReleaseBrief):
            release_brief = result.value
            print("Agent:")
            print(f"Feature: {release_brief.feature}")
            print(f"Benefit: {release_brief.benefit}")
            print(f"Launch date: {release_brief.launch_date}")


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
