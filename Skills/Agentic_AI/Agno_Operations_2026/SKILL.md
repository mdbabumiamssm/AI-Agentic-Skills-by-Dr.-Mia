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
name: agno-operations-2026
description: Build and operate Agno agents, teams, and workflows with state, tools, memory, knowledge, approvals, and production runtime patterns. Use when you want a Python-first agent runtime that tightly connects build-time orchestration to production serving.
measurable_outcome: Select an Agno agent, team, or workflow pattern plus a deployment and governance stance for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Agno Operations (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Agno agent, team, workflow, and production-runtime surfaces.
2. Choose the right abstraction first: single agent, coordinated team, or explicit workflow.
3. Define memory, knowledge, storage, and approval requirements before adding many tools.
4. Decide whether the workload stays local or is served through AgentOS and a monitored runtime.
5. Add a smoke test that covers one stateful turn, one tool path, and one approval or tracing path.

## When to Use

- You want agents, teams, and workflows on one Python runtime.
- Stateful sessions, agentic RAG, or approval workflows matter from day one.
- A local prototype needs a credible path to a served production API.
- MCP-backed tools or shared knowledge layers should be part of the design.

## Output Requirements

- Return the selected Agno abstraction and why.
- State the storage or session model.
- Include one tool or knowledge design note.
- Include one production or governance guardrail.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
