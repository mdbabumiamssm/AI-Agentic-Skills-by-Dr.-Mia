<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Running a simple agent with Anthropic

This sample demonstrates how to create and run a basic agent with Anthropic Claude models.

## What this sample demonstrates

- Creating an AI agent with Anthropic Claude
- Running a simple agent with instructions
- Managing agent lifecycle

## Prerequisites

Before you begin, ensure you have the following prerequisites:

- .NET 8.0 SDK or later
- Anthropic API key configured

**Note**: This sample uses Anthropic Claude models. For more information, see [Anthropic documentation](https://docs.anthropic.com/).

Set the following environment variables:

```powershell
$env:ANTHROPIC_API_KEY="your-anthropic-api-key"  # Replace with your Anthropic API key
$env:ANTHROPIC_MODEL="your-anthropic-model"  # Replace with your Anthropic model
```

## Run the sample

Navigate to the AgentWithAnthropic sample directory and run:

```powershell
cd dotnet\samples\GettingStarted\AgentWithAnthropic
dotnet run --project .\Agent_Anthropic_Step01_Running
```

## Expected behavior

The sample will:

1. Create an agent with Anthropic Claude
2. Run the agent with a simple prompt
3. Display the agent's response



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->