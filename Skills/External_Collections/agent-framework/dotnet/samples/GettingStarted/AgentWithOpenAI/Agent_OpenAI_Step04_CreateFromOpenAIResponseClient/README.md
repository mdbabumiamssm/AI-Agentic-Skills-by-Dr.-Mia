<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Creating an Agent from an OpenAIResponseClient

This sample demonstrates how to create an AI agent directly from an `OpenAI.Responses.OpenAIResponseClient` instance using the `OpenAIResponseClientAgent` class.

## What This Sample Shows

- **Direct OpenAIResponseClient Creation**: Shows how to create an `OpenAI.Responses.OpenAIResponseClient` from `OpenAI.OpenAIClient` and then use it to instantiate an agent
- **OpenAIResponseClientAgent**: Demonstrates using the OpenAI SDK primitives instead of the ones from Microsoft.Extensions.AI and Microsoft.Agents.AI abstractions
- **Full Agent Capabilities**: Shows both regular and streaming invocation of the agent

## Running the Sample

1. Set the required environment variables:
   ```bash
   set OPENAI_API_KEY=your_api_key_here
   set OPENAI_MODEL=gpt-4o-mini
   ```

2. Run the sample:
   ```bash
   dotnet run
   ```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->