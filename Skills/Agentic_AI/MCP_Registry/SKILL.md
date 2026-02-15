---
name: MCP_Registry_Agent
description: A meta-agent that manages Model Context Protocol (MCP) servers and connects them to other agents.
version: 1.0.0
category: Agentic_AI
author: AI Agentic Skills Team
input:
  command: "list", "connect <server_name>", "start <server_name>"
output:
  status: Connection status or server list.
---

# MCP Registry Agent

The **Model Context Protocol (MCP)** is the new standard for connecting AI models to data sources (GitHub, Slack, PostgreSQL). This agent acts as a "Service Discovery" layer for the Agentic OS.

## Capabilities
1.  **Server Management:** Start/Stop local MCP servers (stdio).
2.  **Tool Discovery:** List available tools provided by connected MCP servers.
3.  **Universal Connector:** Allows Claude, Gemini, or DeepSeek agents to use the same toolset.

## Usage
```python
from Skills.Agentic_AI.MCP_Registry.registry import MCPRegistry

registry = MCPRegistry()
registry.start_server("github-mcp")
tools = registry.get_tools("github-mcp")
```
