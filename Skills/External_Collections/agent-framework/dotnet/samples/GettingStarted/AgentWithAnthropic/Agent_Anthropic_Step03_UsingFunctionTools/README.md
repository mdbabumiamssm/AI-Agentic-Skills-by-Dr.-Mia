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

# Using Function Tools with Anthropic agents

This sample demonstrates how to use function tools with Anthropic Claude agents, allowing agents to call custom functions to retrieve information.

## What this sample demonstrates

- Creating function tools using AIFunctionFactory
- Passing function tools to an Anthropic Claude agent
- Running agents with function tools (text output)
- Running agents with function tools (streaming output)
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
dotnet run --project .\Agent_Anthropic_Step03_UsingFunctionTools
```

## Expected behavior

The sample will:

1. Create an agent named "WeatherAssistant" with a GetWeather function tool
2. Run the agent with a text prompt asking about weather
3. The agent will invoke the GetWeather function tool to retrieve weather information
4. Run the agent again with streaming to display the response as it's generated
5. Clean up resources by deleting the agent



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->