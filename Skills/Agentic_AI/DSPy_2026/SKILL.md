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
name: dspy-2026
description: Build and optimize modular LM programs with DSPy. Use when evaluation, compilation, prompt optimization, and metric-driven improvement matter more than adopting a full agent hosting runtime.
measurable_outcome: Select a DSPy module graph, optimization strategy, and evaluation metric for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# DSPy (2026)

## Workflow

1. Read `references/sources.md` and confirm the current DSPy guidance for modules, optimizers, MCP usage, observability, and deployment.
2. Decide whether the need is prompt optimization, modular LM programming, RAG improvement, or an agent loop that benefits from compilation and evaluation.
3. Define the signature and evaluation metric before selecting an optimizer.
4. Keep runtime orchestration separate from optimization concerns: DSPy can power an agent layer without becoming the whole hosting platform.
5. Add a smoke test that covers one compiled program run, one evaluation metric, and one optimizer or fallback decision.

## When to Use

- Prompt strings are getting brittle and need metric-driven optimization.
- The system should stay modular and portable across models.
- You want evaluation-first LM engineering rather than ad hoc prompt iteration.
- An agent or RAG pipeline needs compilation, signatures, and measurable improvement loops.

## Output Requirements

- Return the selected DSPy module pattern and why.
- State the optimization target and evaluation metric.
- Include one note on how DSPy coexists with the runtime framework.
- Include one debugging, observability, or rollback note.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
