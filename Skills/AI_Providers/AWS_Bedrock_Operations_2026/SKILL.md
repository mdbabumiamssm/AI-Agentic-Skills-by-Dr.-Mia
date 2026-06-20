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

- Review the Strands Agents and Amazon Bedrock AgentCore reference architecture for MCP and A2A integration, browser automation and voice tools, explicit identity boundaries, end-to-end observability, deployment topology, and production hardening requirements.
- Validate Strands Agents plus Amazon Bedrock AgentCore reference architectures for multimodal or voice-enabled agents against explicit identity and session-state ownership, MCP and A2A trust boundaries, isolated browser automation, end-to-end observability, documented deployment topology, and least-privilege controls, using `aws-samples/sample-strands-agent-with-agentcore` as the implementation reference.
- Use the 2026-05-04 `aws-samples/sample-strands-agent-with-agentcore` TypeScript reference architecture as an implementation example for AWS-native agentic chatbots with Strands Agents and Amazon Bedrock AgentCore; map Strands agent orchestration, Bedrock AgentCore deployment/runtime, MCP integration, agent-to-agent (A2A) protocol use, browser automation and voice assistant patterns, identity and session boundaries, least-privilege IAM, TypeScript deployment topology, production observability checks, and production guardrails before adapting it.
- Adapt the Strands Agents on Bedrock AgentCore sample by defining MCP and A2A trust boundaries, workload identity and least-privilege IAM, isolated browser and voice tool runtimes, observable agent and tool paths, explicit session-state ownership and retention, deployment and rollback topology, human approval gates for consequential actions, and production hardening before deployment.

## Strands Agents + AgentCore Reference Architecture

- Treat the sample as a TypeScript chatbot reference: keep chat or channel entrypoints, Strands agent orchestration, AgentCore hosting/runtime binding, Bedrock access, and tool adapters as explicit boundaries.
- Review MCP tool contracts and A2A handoff paths before adding tools, multi-agent routing, browser automation, or voice assistant channels.
- For browser automation or voice assistant tools, document the invocation path, permission boundary, and channel-specific failure handling before enabling them in an agent flow.
- Keep browser automation and voice assistant patterns behind explicit permission, network, and channel-handling boundaries.
- Define least-privilege IAM boundaries for each agent, tool, browser automation path, and channel surface before deployment.
- Document the deployment topology across chatbot ingress, Strands agent code, AgentCore runtime, Bedrock model access, MCP/browser/voice integrations, logs, metrics, and rollback units.
- Use AgentCore when the workload should follow the reference architecture's hosted agent runtime pattern; use custom orchestration when existing routing, deployment, or control-plane requirements make that hosted topology a poor fit, and document the decision before implementation.
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
