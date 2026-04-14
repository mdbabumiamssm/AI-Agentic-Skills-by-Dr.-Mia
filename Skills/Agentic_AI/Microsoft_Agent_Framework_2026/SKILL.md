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
name: microsoft-agent-framework-2026
description: Build and operate Microsoft Agent Framework workflows across Python and .NET. Use when enterprise hosting, workflow graphs, or migration from older Microsoft agent stacks matters.
measurable_outcome: Select a Microsoft Agent Framework workflow pattern, hosting direction, and migration stance for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Microsoft Agent Framework (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Agent Framework guidance before starting implementation.
2. Decide whether the workload is a simple agent, a workflow, or a multi-agent composition.
3. Record the hosting target and the language boundary early: Python, .NET, or mixed estate.
4. If an existing stack uses AutoGen or Semantic Kernel, write the migration posture explicitly instead of blending abstractions.
5. Add a smoke test that covers conversation state, one tool path, and one workflow transition.

## When to Use

- The team spans Python and .NET.
- Enterprise hosting and lifecycle management matter.
- You need workflow graphs and persistent conversational state.
- You are consolidating prior AutoGen or Semantic Kernel work.

## Output Requirements

- Return the chosen Agent Framework pattern and language path.
- Include one migration note for older Microsoft agent stacks when relevant.
- Include one hosting or deployment assumption.
- Include one workflow-level smoke test.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
