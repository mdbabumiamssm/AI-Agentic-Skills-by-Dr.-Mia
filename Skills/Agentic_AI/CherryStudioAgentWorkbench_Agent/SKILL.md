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
name: 'cherry-studio-agent-workbench'
description: 'Operate Cherry Studio as a multi-provider AI workbench for assistants, autonomous agents, coding-agent integrations, local model routing, credentials, and guardrails.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Cherry Studio Agent Workbench

## Overview

Use this skill to configure, operate, and troubleshoot Cherry Studio as a desktop-style AI productivity workbench for multi-provider chat, assistants, autonomous agents, and coding-agent integrations. It helps agents make practical decisions about provider routing, credentials, local model use, assistant libraries, skill workflows, and governance controls without treating any single model runtime as the whole system.

## When to Use This Skill

- Setting up or validating Cherry Studio as a multi-provider AI assistant and agent workbench.
- Configuring assistants, autonomous-agent workflows, skills, or prompt libraries inside Cherry Studio.
- Routing tasks across frontier LLM providers, local models, coding agents, or CLI-backed agent tools.
- Troubleshooting provider credentials, model availability, agent tool access, workspace behavior, or integration issues.
- Designing guardrails for team usage, sensitive data handling, auditability, and cost-aware model selection.
- Comparing Cherry Studio workbench behavior with direct provider APIs or standalone coding-agent tools.

## Core Capabilities

1. Provider and model routing: Map user tasks to suitable hosted or local models, then document provider, model, context, cost, privacy, and latency tradeoffs.
2. Assistant workspace design: Create or refine assistants with clear roles, instructions, tools, knowledge sources, and escalation boundaries.
3. Autonomous-agent operation: Plan agent workflows with explicit goals, allowed actions, checkpointing, human review points, and rollback expectations.
4. Coding-agent integration: Coordinate Cherry Studio usage with coding-oriented tools such as Codex, Claude Code, OpenCode, OpenCLI, or similar integrations when present in the user environment.
5. Credential and secret handling: Verify that API keys, local endpoints, and environment variables are configured without exposing secrets in prompts, logs, screenshots, or shared artifacts.
6. Local model workflows: Route appropriate tasks to local inference endpoints when privacy, offline use, or provider independence matters, while recording known limitations.
7. Governance guardrails: Define policies for sensitive data, model selection, tool permissions, workspace sharing, logging, and review of autonomous actions.
8. Troubleshooting and validation: Use a minimal reproducible prompt, provider status checks, configuration inspection, and a short end-to-end test before declaring the workflow ready.

## Inputs / Outputs

Inputs consumed by this skill:

- User goal, task type, risk level, and expected output format.
- Cherry Studio version or repository state when available.
- Provider list, model preferences, local model endpoints, and credential configuration status.
- Assistant definitions, system prompts, skill instructions, tool permissions, or workspace settings.
- Error messages, logs, screenshots, or reproduction steps for troubleshooting.

Outputs produced by this skill:

- A Cherry Studio setup or operating plan with provider routing and assistant/workflow recommendations.
- A validated assistant or autonomous-agent configuration checklist.
- A troubleshooting diagnosis with concrete remediation steps.
- A governance note covering secrets, sensitive data, tool permissions, auditability, and human review points.
- A brief validation result showing the prompt used, selected model/provider, observed output, and remaining risks.

## References

- CherryHQ/cherry-studio GitHub repository: https://github.com/CherryHQ/cherry-studio
