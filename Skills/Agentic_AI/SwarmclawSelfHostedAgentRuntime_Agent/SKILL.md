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



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'swarmclaw-self-hosted-agent-runtime'
description: 'Deploy and operate SwarmClaw as a self-hosted multi-agent AI runtime with MCP servers, memory, skills, schedules, and multi-provider LLM routing.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# SwarmClaw Self-Hosted Agent Runtime

## Overview

This skill guides setup, configuration, and operation of SwarmClaw, a self-hosted multi-agent AI runtime with MCP server support, memory, skills, schedules, and integrations across many LLM providers. Use it to turn a SwarmClaw repository or deployment target into a runnable agent runtime with clear provider boundaries, orchestration patterns, and operational checks.

## When to Use This Skill

- Deploying or evaluating the `swarmclawai/swarmclaw` runtime.
- Configuring self-hosted multi-agent orchestration with MCP servers.
- Wiring LLM providers such as Claude, GPT, Gemini, OpenRouter, or Ollama into one runtime.
- Designing agent schedules, memory usage, skill execution, and autonomous task boundaries.
- Reviewing security, secrets, network access, and tool permission boundaries for a self-hosted agent platform.
- Debugging runtime startup, provider selection, MCP connectivity, or agent execution flow.

## Core Capabilities

1. **Runtime orientation** - Inspect the SwarmClaw repository, identify service boundaries, required environment variables, startup commands, and deployment assumptions before making changes.
2. **Provider configuration** - Map requested models to supported providers, verify required API keys or local endpoints, and keep provider-specific secrets out of source-controlled files.
3. **MCP integration** - Configure MCP servers, validate tool discovery, and confirm that each agent only receives the tools needed for its role.
4. **Agent orchestration design** - Define agents, skills, memory policies, schedules, and task handoffs with explicit ownership and stop conditions.
5. **Self-hosted deployment checks** - Validate local or server deployment prerequisites, process management, logs, ports, health checks, and persistent storage.
6. **Security boundary review** - Check secret handling, network exposure, file-system access, tool permissions, scheduler privileges, and auditability before enabling autonomous execution.
7. **Operational troubleshooting** - Use logs, configuration diffs, provider test calls, MCP probes, and minimal reproduction tasks to isolate runtime failures.

## Inputs / Outputs

**Inputs**

- SwarmClaw repository path or deployment target.
- Intended runtime mode: local development, private server, containerized service, or production deployment.
- LLM provider choices, model names, API keys, base URLs, and local inference endpoints.
- MCP server definitions, allowed tools, memory/storage requirements, and scheduling needs.
- Agent roles, task objectives, autonomy limits, and security constraints.

**Outputs**

- A validated SwarmClaw setup or implementation plan.
- Provider and MCP configuration guidance tailored to the target environment.
- Agent orchestration design covering roles, memory, skills, schedules, and handoff rules.
- Deployment checklist with startup, health-check, logging, and persistence steps.
- Security review notes for secrets, permissions, network exposure, and autonomous execution controls.
- Troubleshooting summary with observed failures, likely causes, and next commands or patches.

## References

- SwarmClaw GitHub repository: https://github.com/swarmclawai/swarmclaw
