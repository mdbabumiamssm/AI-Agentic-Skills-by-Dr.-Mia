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
name: langgraph-production-2026
description: Build durable LangGraph-based agents with checkpoints, interrupts, human review, and observability. Use when agent state, resumability, or long-running orchestration matter more than lightweight demos.
measurable_outcome: Design or update a production LangGraph workflow with persistence, interrupt points, and a minimal smoke-test plan within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# LangGraph Production (2026)

## Workflow

1. Read `references/sources.md` and confirm the current LangGraph and LangSmith primitives before coding.
2. Decide whether the workload should be a graph, a deep-agent coordinator, or a simpler provider-native agent.
3. Add persistence and checkpoints first, then define interrupt and resume boundaries for risky tool calls.
4. Add human-in-the-loop review where write actions, execution, or data mutation can occur.
5. Define tracing and evaluation hooks before expanding the workflow surface.

## When to Use

- Durable execution and resumability matter.
- The system needs explicit graph state rather than hidden prompt state.
- Human review must pause and resume tool execution safely.
- Multi-step workflows may run longer than a single synchronous request.

## Output Requirements

- Return the selected graph pattern with one-line rationale.
- State where checkpoints are stored and where interrupts occur.
- Include one LangSmith tracing path or equivalent observability note.
- Include a minimal smoke-test plan covering resume, tool failure, and approval flow.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
