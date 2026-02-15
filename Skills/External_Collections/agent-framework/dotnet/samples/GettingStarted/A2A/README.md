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

# Agent-to-Agent (A2A) Samples

These samples demonstrate how to work with Agent-to-Agent (A2A) specific features in the Agent Framework.

For other samples that demonstrate how to use AIAgent instances,
see the [Getting Started With Agents](../Agents/README.md) samples.

## Prerequisites

See the README.md for each sample for the prerequisites for that sample.

## Samples

|Sample|Description|
|---|---|
|[A2A Agent As Function Tools](./A2AAgent_AsFunctionTools/)|This sample demonstrates how to represent an A2A agent as a set of function tools, where each function tool corresponds to a skill of the A2A agent, and register these function tools with another AI agent so it can leverage the A2A agent's skills.|
|[A2A Agent Polling For Task Completion](./A2AAgent_PollingForTaskCompletion/)|This sample demonstrates how to poll for long-running task completion using continuation tokens with an A2A agent.|

## Running the samples from the console

To run the samples, navigate to the desired sample directory, e.g.

```powershell
cd A2AAgent_AsFunctionTools
```

Set the required environment variables as documented in the sample readme.
If the variables are not set, you will be prompted for the values when running the samples.
Execute the following command to build the sample:

```powershell
dotnet build
```

Execute the following command to run the sample:

```powershell
dotnet run --no-build
```

Or just build and run in one step:

```powershell
dotnet run
```

## Running the samples from Visual Studio

Open the solution in Visual Studio and set the desired sample project as the startup project. Then, run the project using the built-in debugger or by pressing `F5`.

You will be prompted for any required environment variables if they are not already set.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->