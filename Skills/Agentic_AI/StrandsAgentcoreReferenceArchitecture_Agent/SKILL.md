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
name: 'strands-agentcore-reference-architecture'
description: 'Implement AWS Strands Agents with Amazon Bedrock AgentCore reference architectures for chat, voice, A2A/MCP, browser automation, observability, deployment, and security.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Strands AgentCore Reference Architecture

## Overview

Use this skill to design, implement, or review agentic AI chatbot and assistant systems based on AWS Strands Agents with Amazon Bedrock AgentCore. It focuses on practical reference-architecture decisions: agent runtime shape, tool interoperability, chat or voice interfaces, deployment boundaries, observability, and security controls.

## When to Use This Skill

- Implementing a Strands Agents application on Amazon Bedrock AgentCore.
- Adapting the `aws-samples/sample-strands-agent-with-agentcore` architecture for a production or prototype assistant.
- Designing agent-to-agent interoperability with A2A protocol concepts.
- Exposing or consuming tools through MCP servers or MCP-compatible clients.
- Adding browser automation, web task execution, or interactive tool use to an agent workflow.
- Building chat or voice assistant surfaces backed by Bedrock-hosted agent runtime components.
- Reviewing IAM, network, credential, audit, and observability boundaries for an AWS agent runtime.
- Translating a TypeScript sample architecture into deployable infrastructure, application modules, or CI/CD tasks.

## Core Capabilities

1. **Reference architecture mapping** - Identify the user-facing surface, orchestration layer, Strands agent runtime, Bedrock AgentCore integration, tool services, state stores, and deployment units before coding.

2. **AWS sample reference implementation** - Use `aws-samples/sample-strands-agent-with-agentcore` as the TypeScript baseline for Strands Agents plus Bedrock AgentCore patterns across A2A/MCP integration, browser automation, voice assistant paths, deployment boundaries, IAM/security controls, and observability handoff.

3. **TypeScript implementation planning** - Preserve the repository's TypeScript deployment shape, package boundaries, environment configuration, and build scripts when extending the sample.

4. **AgentCore runtime integration** - Connect agent handlers, model invocation paths, memory or session state, and runtime configuration to Amazon Bedrock AgentCore patterns.

5. **A2A interoperability design** - Define agent cards, capability metadata, request routing, identity expectations, and handoff behavior when agents need to collaborate.

6. **MCP interoperability design** - Treat MCP servers as explicit tool boundaries with typed inputs, scoped permissions, error handling, and auditable tool-call traces.

7. **Browser automation support** - Isolate browser automation tools from core agent logic, constrain browsing permissions, capture durable evidence, and handle timeouts or nondeterministic page behavior.

8. **Chat and voice surface guidance** - Separate transport concerns from agent reasoning so web chat, streaming chat, and voice assistant flows can share orchestration logic.

9. **Observability and evaluation** - Instrument requests, model calls, tool calls, traces, latency, failures, session IDs, and human feedback without logging secrets or sensitive payloads.

10. **IAM and security boundaries** - Apply least-privilege IAM, scoped secrets, per-tool authorization, network egress controls, data-retention limits, and clear trust boundaries between agents, tools, and users.

11. **Deployment readiness review** - Check environment variables, infrastructure outputs, rollback approach, alarms, dependency versions, and operational runbooks before promoting changes.

## Inputs / Outputs

**Inputs**

- User goal, target assistant surface, and expected agent behaviors.
- Existing repository files, deployment manifests, package scripts, and environment configuration.
- Required AWS account, region, Bedrock model, AgentCore runtime, IAM, logging, and network constraints.
- Tooling requirements for A2A, MCP, browser automation, external APIs, chat, or voice workflows.
- Security, privacy, compliance, and observability requirements for the intended deployment.

**Outputs**

- A concise architecture plan showing components, responsibilities, data flow, and trust boundaries.
- TypeScript implementation guidance or code changes aligned with the repository structure.
- A tool interoperability plan for A2A and MCP capabilities, including error and permission handling.
- A deployment checklist covering environment variables, IAM permissions, observability, rollback, and validation.
- Security review notes that identify risks, mitigations, and assumptions needing user or platform-owner confirmation.

## References

- Source reference architecture: https://github.com/aws-samples/sample-strands-agent-with-agentcore
