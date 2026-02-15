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

# A2A Agent as Function Tools

This sample demonstrates how to represent an A2A agent as a set of function tools, where each function tool corresponds to a skill of the A2A agent, 
and register these function tools with another AI agent so it can leverage the A2A agent's skills.

# Prerequisites

Before you begin, ensure you have the following prerequisites:

- .NET 10 SDK or later
- Access to the A2A agent host service

**Note**: These samples need to be run against a valid A2A server. If no A2A server is available, they can be run against the echo-agent that can be 
spun up locally by following the guidelines at: https://github.com/a2aproject/a2a-dotnet/blob/main/samples/AgentServer/README.md

Set the following environment variables:

```powershell
$env:A2A_AGENT_HOST="https://your-a2a-agent-host" # Replace with your A2A agent host endpoint
$env:AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com/" # Replace with your Azure OpenAI resource endpoint
$env:AZURE_OPENAI_DEPLOYMENT_NAME="gpt-4o-mini"  # Optional, defaults to gpt-4o-mini
```

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->