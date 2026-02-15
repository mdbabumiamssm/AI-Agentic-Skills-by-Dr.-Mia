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

# Multi-Agent Sample

This sample demonstrates how to use the Durable Extension for Agent Framework to create an Azure Functions app that hosts multiple AI agents and provides direct HTTP API access for interactive conversations with each agent.

## Key Concepts Demonstrated

- Using the Microsoft Agent Framework to define multiple AI agents with unique names and instructions.
- Registering multiple agents with the Function app and running them using HTTP.
- Conversation management (via thread IDs) for isolated interactions per agent.
- Two different methods for registering agents: list-based initialization and incremental addition.

## Prerequisites

Complete the common environment preparation steps described in `../README.md`, including installing Azure Functions Core Tools, starting Azurite, configuring Azure OpenAI settings, and installing this sample's requirements.

## Running the Sample

With the environment setup and function app running, you can test the sample by sending HTTP requests to the different agent endpoints.

You can use the `demo.http` file to send messages to the agents, or a command line tool like `curl` as shown below:

> **Note:** Each endpoint waits for the agent response by default. To receive an immediate HTTP 202 instead, set the `x-ms-wait-for-response` header or include `"wait_for_response": false` in the request body.

### Test the Weather Agent

Bash (Linux/macOS/WSL):
Weather agent request:

```bash
curl -X POST http://localhost:7071/api/agents/WeatherAgent/run \
    -H "Content-Type: application/json" \
    -d '{"message": "What is the weather in Seattle?"}'
```

Expected HTTP 202 payload:

```json
{
  "status": "accepted",
  "response": "Agent request accepted",
  "message": "What is the weather in Seattle?",
  "thread_id": "<guid>",
  "correlation_id": "<guid>"
}
```

Math agent request:

```bash
curl -X POST http://localhost:7071/api/agents/MathAgent/run \
    -H "Content-Type: application/json" \
    -d '{"message": "Calculate a 20% tip on a $50 bill"}'
```

Expected HTTP 202 payload:

```json
{
  "status": "accepted",
  "response": "Agent request accepted",
  "message": "Calculate a 20% tip on a $50 bill",
  "thread_id": "<guid>",
  "correlation_id": "<guid>"
}
```

Health check (optional):

```bash
curl http://localhost:7071/api/health
```

Expected response:

```json
{
  "status": "healthy",
  "agents": [
    {"name": "WeatherAgent", "type": "ChatAgent"},
    {"name": "MathAgent", "type": "ChatAgent"}
  ],
  "agent_count": 2
}
```

## Code Structure

The sample demonstrates two ways to register multiple agents:

### Option 1: Pass list of agents during initialization
```python
app = AgentFunctionApp(agents=[weather_agent, math_agent])
```

### Option 2: Add agents incrementally (commented in sample)
```python
app = AgentFunctionApp()
app.add_agent(weather_agent)
app.add_agent(math_agent)
```

Each agent automatically gets:
- `POST /api/agents/{agent_name}/run` - Send messages to the agent



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->