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
name: openai-agents-sdk-2026
description: Build production agents on OpenAI using the Responses API, Agents SDK, guardrails, built-in tools, tracing, and evals. Use when OpenAI is the primary runtime and you need a modern agent stack rather than legacy Assistants-era patterns.
measurable_outcome: Select an OpenAI agent runtime pattern, tool surface, guardrail path, and validation plan for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# OpenAI Agents SDK (2026)

## Workflow

1. Read `references/sources.md` and confirm the current OpenAI agent platform surfaces before implementation.
2. Choose the runtime path deliberately: Responses-only, full Agents SDK orchestration, ChatKit UI integration, or Agent Builder/hosted workflows.
3. Define the tool surface first: built-in tools, MCP-backed tools, custom function tools, or a mixed design.
4. Add guardrails, tracing, and evals before broadening agent autonomy or adding side effects.
5. Run a smoke test that covers one streamed turn, one tool call, one guarded path, and one evaluation check.

## When to Use

- OpenAI is the primary provider and first-party agent tooling matters.
- You need built-in tools such as web search, file search, computer use, or MCP-backed integrations.
- The workflow needs first-party tracing, evals, background execution, or UI embedding.
- You are migrating from older Assistants-era patterns and want the current OpenAI stack.

## Output Requirements

- Return the chosen OpenAI runtime pattern and why it fits.
- State which tool surface is primary and whether MCP is involved.
- Include one guardrail or approval boundary.
- Include one tracing or evaluation note for rollout safety.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
