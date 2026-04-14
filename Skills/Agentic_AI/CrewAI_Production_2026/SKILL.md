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
name: crewai-production-2026
description: Build production-oriented CrewAI systems using agents, crews, and flows with memory, knowledge, guardrails, and tracing. Use when role-based multi-agent automation is central and you want a dedicated Python framework for crews plus event-driven flow control.
measurable_outcome: Select a CrewAI crew or flow pattern, observability path, and rollout guardrails for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# CrewAI Production (2026)

## Workflow

1. Read `references/sources.md` and confirm the current CrewAI concepts for agents, crews, flows, and tracing.
2. Choose whether the workload should be a crew, a flow, or an intentional combination of both.
3. Define memory, knowledge, and structured output expectations before scaling the agent count.
4. Add tracing and evaluation visibility before rollout, especially if multiple agents can hand work to each other.
5. Run a smoke test that covers one multi-agent path, one flow transition, and one failure or retry condition.

## When to Use

- Role-based multi-agent collaboration is the core operating model.
- You want clear separation between autonomous crews and explicit event-driven flows.
- The team prefers a Python-native framework with first-party tracing and quickstarts.
- The workload benefits from memory, knowledge, and structured outputs without building the orchestration layer from scratch.

## Output Requirements

- Return the selected CrewAI pattern and why.
- State whether crews, flows, or both are used.
- Include one tracing or observability note.
- Include one rollout guardrail for retries, approvals, or telemetry.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
