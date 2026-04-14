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
name: llamaindex-workflows-2026
description: Build event-driven LlamaIndex workflows for document and data agents with serving, human-in-the-loop checkpoints, and workflow observability. Use when retrieval-heavy or document-centric systems need explicit workflow control.
measurable_outcome: Select a LlamaIndex workflow pattern, data boundary, and serving path for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# LlamaIndex Workflows (2026)

## Workflow

1. Read `references/sources.md` and confirm the current workflow, serving, and observability guidance.
2. Decide whether the workload is a retrieval-heavy data agent, an agent workflow, or a document-processing pipeline with agent steps.
3. Define the workflow events, shared state, and data boundaries before tuning prompts or retrieval.
4. Choose the serving path early: local Python process, workflow server, or broader LlamaIndex deployment.
5. Add a smoke test that covers one event transition, one retrieval or tool path, and one visibility or retry signal.

## When to Use

- The agent is tightly coupled to data, indexes, documents, or retrieval workflows.
- Event-driven control is needed instead of a hidden prompt loop.
- You want a server path for workflows without building the orchestration layer from scratch.
- Human-in-the-loop or durable execution should be possible around document-aware tasks.

## Output Requirements

- Return the selected workflow shape and why.
- State the primary data or index boundary.
- Include one serving or deployment note.
- Include one observability or workflow-debugging note.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
