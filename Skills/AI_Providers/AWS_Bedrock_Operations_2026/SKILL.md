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

---
name: aws-bedrock-operations-2026
description: Build and run production workloads on Amazon Bedrock with current model access, IAM controls, and reliability guardrails. Use when implementing Bedrock inference pipelines, agent workflows, or model upgrades.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# AWS Bedrock Operations (2026)

## Core Capabilities

- When adapting the 2026-05-04 `aws-samples/sample-strands-agent-with-agentcore` TypeScript reference architecture, document the Amazon Bedrock AgentCore hosting model, MCP tool integration, browser automation and voice-assistant surfaces, A2A/agent-to-agent assumptions, IAM boundaries for agents, tools, and channels, and deployment review checkpoints before production use.
- Use `aws-samples/sample-strands-agent-with-agentcore` as a 2026-05-04 TypeScript reference architecture for Strands Agents plus Amazon Bedrock AgentCore agentic chatbots: map MCP and A2A integration points, isolate browser automation and voice assistant patterns, define IAM boundaries around agents and tools, and require production guardrails before deployment.
- Use the 2026-05-04 `aws-samples/sample-strands-agent-with-agentcore` TypeScript reference architecture as an implementation example for agentic chatbots with Strands Agents and Amazon Bedrock AgentCore, covering A2A protocol flows, MCP tool integration, browser automation, voice assistant channels, IAM least-privilege boundaries, deployment review, and production guardrail checks before adapting it.
- Use `aws-samples/sample-strands-agent-with-agentcore` as a reference architecture for Strands Agents plus Amazon Bedrock AgentCore agentic chatbots: verify MCP/A2A interoperability contracts, browser and voice assistant components, IAM and networking guardrails, deployment topology, and the gap between sample scaffolding and production-ready infrastructure before adoption.
- Use the `aws-samples/sample-strands-agent-with-agentcore` TypeScript sample as source-backed architecture input for Strands Agents plus Amazon Bedrock AgentCore: capture AgentCore deployment boundaries, MCP/A2A integration paths, browser and voice assistant channel patterns, IAM isolation, observability expectations, and production promotion checklist items before adapting it.
- Use `aws-samples/sample-strands-agent-with-agentcore` as a reference pattern for Strands Agents with Amazon Bedrock AgentCore chatbots: map chatbot architecture, MCP and browser automation hooks, A2A-style orchestration, voice assistant surfaces, IAM boundaries, and deployment hygiene before adapting the TypeScript sample to a target environment.
- Use the `aws-samples/sample-strands-agent-with-agentcore` sample as architecture guidance for Strands Agents plus Amazon Bedrock AgentCore designs: inventory MCP integrations, A2A-style agent communication paths, browser automation and voice-assistant surfaces, the TypeScript deployment shape, and IAM boundaries; treat sample code as reference architecture rather than production-ready code until reviewed for the target environment.
- Treat `aws-samples/sample-strands-agent-with-agentcore` as a Strands Agents with Amazon Bedrock AgentCore reference architecture for agentic chatbot prototypes: evaluate AgentCore hosting topology, MCP/tool integration, browser automation or voice-assistant channels, and IAM boundaries before using the sample as production guidance.
- Use `aws-samples/sample-strands-agent-with-agentcore` (published 2026-05-04) as a TypeScript reference architecture when designing Strands Agents with Amazon Bedrock AgentCore: review MCP tool integration, browser automation permissions, voice and chatbot channel patterns, A2A handoff considerations, IAM role boundaries, and deployment validation steps before adapting the sample for production agent workflows.
- Use Strands Agents plus Amazon Bedrock AgentCore reference architectures as an implementation pattern for managed Bedrock agentic chatbots, including AgentCore runtime boundaries and deployment topology, A2A and MCP integration, browser automation and voice/chat channel surfaces, IAM isolation and least privilege, observability hooks, production chatbot guardrails, evaluation gates, and deployment caveats such as service limits and whether the TypeScript sample shape fits the target workflow.
- Use the 2026-05-04 `aws-samples/sample-strands-agent-with-agentcore` GitHub reference architecture as a Strands Agents plus Amazon Bedrock AgentCore chatbot pattern: map A2A and MCP integrations explicitly, keep browser automation behind clear permission and network boundaries, separate voice-assistant channel handling from chat orchestration, isolate agent and tool IAM roles, document runtime and deployment topology, and treat the TypeScript sample repo as reference scaffolding until reviewed against the target environment.
- Use the `aws-samples/sample-strands-agent-with-agentcore` architecture as a reference point for MCP tool integration, agent-to-agent patterns, browser automation and voice assistant components, IAM boundary reviews, and deployment topology planning; treat its TypeScript sample code as non-production scaffolding until service limits, security controls, observability, tests, and operational runbooks are validated for the target environment.
- For Strands Agents with Amazon Bedrock AgentCore chatbot designs based on `aws-samples/sample-strands-agent-with-agentcore`, document MCP tool contracts, A2A handoff paths, browser or voice-assistant surfaces, per-agent and tool IAM boundaries, observability signals, and sample-app caveats before promoting the architecture beyond reference use.
- Evaluate AgentCore chatbot designs against the Strands Agents, Amazon Bedrock, A2A protocol, MCP, browser automation, and voice assistant patterns demonstrated by `aws-samples/sample-strands-agent-with-agentcore`.
- Use `aws-samples/sample-strands-agent-with-agentcore` as source-backed guidance for Strands Agents plus Amazon Bedrock AgentCore chatbot and voice-agent deployments: define MCP integration points, keep browser automation behind explicit permission and network boundaries, separate agent and tool IAM roles, instrument agent handoffs and tool/channel flows, and require guardrails, tests, rollback criteria, and operational review before production use.
- For multi-agent chatbots using Strands Agents plus Amazon Bedrock AgentCore, use `aws-samples/sample-strands-agent-with-agentcore` as reference-architecture input when choosing the AgentCore runtime shape, defining per-agent and tool IAM boundaries, integrating MCP tools, isolating browser or voice assistant components, mapping deployment topology, and setting production guardrails before adapting the TypeScript sample.

## Strands Agents + AgentCore Reference Architecture

- Start from Strands Agents on Amazon Bedrock when the chatbot needs managed Bedrock agent workflows with explicit AgentCore runtime, deployment topology, and rollback boundaries.
- Review MCP tool contracts and A2A handoff paths before adding tools, multi-agent routing, browser automation, or voice assistant channels.
- Keep browser automation and voice assistant patterns behind explicit permission, network, and channel-handling boundaries.
- Define least-privilege IAM boundaries for each agent, tool, browser automation path, and channel surface before deployment.
- Add observability hooks for inference calls, agent handoffs, MCP tool calls, browser automation steps, voice assistant flows, errors, retries, and throttling.
- Production checklist: confirm regional model availability, service limits, guardrails, tests, evaluation gates, security review, operational runbooks, and rollback criteria before treating the TypeScript sample architecture as deployable.

## Workflow

1. Confirm regional model availability, access grants, and service limits from `references/sources.md`.
2. Choose the Bedrock path that matches the workload: Converse API, Agents, Knowledge Bases, Guardrails, or batch evaluation.
3. Implement IAM, secret handling, and network boundaries before application logic.
4. Define retries, throttling backoff, timeout budgets, and observability around every inference call.
5. Validate with a small canary run before broader rollout and record rollback triggers.

## Output Requirements

- State the selected Bedrock model or service path and why.
- State the IAM and network pattern being used.
- Include one operational guardrail for reliability or cost control.

## References

- https://github.com/aws-samples/sample-strands-agent-with-agentcore


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
