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
name: mastra-production-2026
description: Build production-oriented Mastra agents and workflows with memory, MCP, evals, observability, and TypeScript-first application integration. Use when product teams want one TypeScript stack for agent runtime, workflows, and operational tooling.
measurable_outcome: Select a Mastra agent or workflow pattern, runtime boundary, and observability path for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Mastra Production (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Mastra surface for agents, workflows, MCP, evals, and observability.
2. Decide whether the workload should center on an embedded agent, a workflow, or a combined product runtime.
3. Define memory, tool access, and MCP boundaries before wiring many integrations into the app.
4. Choose the runtime path early: local development, Mastra Cloud, or a self-hosted TypeScript service.
5. Add a smoke test that covers one tool path, one workflow transition, and one observability or eval signal.

## When to Use

- The team is TypeScript-first and wants a modern agent runtime inside the application stack.
- Agents, workflows, MCP, evals, and observability should live in one framework.
- The workload needs product-facing agents rather than a research-only orchestration layer.
- You want a credible path from local prototype to deployed runtime without rebuilding the stack.

## Output Requirements

- Return the selected Mastra runtime pattern and why.
- State whether the workload centers on agents, workflows, or both.
- Include one MCP or tool-surface assumption.
- Include one observability or eval note for rollout safety.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
