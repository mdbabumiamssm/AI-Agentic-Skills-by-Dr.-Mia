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

- Evaluate Strands Agents plus Amazon Bedrock AgentCore reference architectures for Bedrock agent workloads, including A2A/MCP interoperability, browser automation boundaries, voice/chatbot patterns, IAM boundaries, production deployment checks, and whether the TypeScript sample shape fits the target workflow.

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
