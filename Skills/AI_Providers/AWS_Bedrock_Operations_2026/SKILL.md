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

- Treat `aws-samples/sample-strands-agent-with-agentcore` as a Strands Agents with Amazon Bedrock AgentCore reference architecture for agentic chatbot prototypes: evaluate AgentCore hosting topology, MCP/tool integration, browser automation or voice-assistant channels, and IAM boundaries before using the sample as production guidance.
- Use `aws-samples/sample-strands-agent-with-agentcore` (published 2026-05-04) as a TypeScript reference architecture when designing Strands Agents with Amazon Bedrock AgentCore: review MCP tool integration, browser automation permissions, voice and chatbot channel patterns, A2A handoff considerations, IAM role boundaries, and deployment validation steps before adapting the sample for production agent workflows.
- Use Strands Agents plus Amazon Bedrock AgentCore reference architectures as an implementation pattern for managed Bedrock agentic chatbots, including AgentCore runtime boundaries and deployment topology, A2A and MCP integration, browser automation and voice/chat channel surfaces, IAM isolation and least privilege, observability hooks, production chatbot guardrails, evaluation gates, and deployment caveats such as service limits and whether the TypeScript sample shape fits the target workflow.
- Use the 2026-05-04 `aws-samples/sample-strands-agent-with-agentcore` GitHub reference architecture as a Strands Agents plus Amazon Bedrock AgentCore chatbot pattern: map A2A and MCP integrations explicitly, keep browser automation behind clear permission and network boundaries, separate voice-assistant channel handling from chat orchestration, isolate agent and tool IAM roles, document runtime and deployment topology, and treat the TypeScript sample repo as reference scaffolding until reviewed against the target environment.
- Use the `aws-samples/sample-strands-agent-with-agentcore` architecture as a reference point for MCP tool integration, agent-to-agent patterns, browser automation and voice assistant components, IAM boundary reviews, and deployment topology planning; treat its TypeScript sample code as non-production scaffolding until service limits, security controls, observability, tests, and operational runbooks are validated for the target environment.
- For Strands Agents with Amazon Bedrock AgentCore chatbot designs based on `aws-samples/sample-strands-agent-with-agentcore`, document MCP tool contracts, A2A handoff paths, browser or voice-assistant surfaces, per-agent and tool IAM boundaries, observability signals, and sample-app caveats before promoting the architecture beyond reference use.
- Evaluate AgentCore chatbot designs against the Strands Agents, Amazon Bedrock, A2A protocol, MCP, browser automation, and voice assistant patterns demonstrated by `aws-samples/sample-strands-agent-with-agentcore`.

## AgentCore Reference Architecture

- Start from Strands Agents on Amazon Bedrock when the chatbot needs managed Bedrock agent workflows with explicit AgentCore runtime and deployment boundaries.
- Review A2A and MCP integration points before adding tools, multi-agent handoffs, browser automation, or voice assistant channels.
- Define IAM boundaries for each agent, tool, browser automation path, and channel surface before deployment.
- Apply guardrails, observability hooks, tests, and operational runbooks before treating the TypeScript sample architecture as production-ready.
- Check regional model availability, service limits, deployment topology, and rollback criteria for the target environment.

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
