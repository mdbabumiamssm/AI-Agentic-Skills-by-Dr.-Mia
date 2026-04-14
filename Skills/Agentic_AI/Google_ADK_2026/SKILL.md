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
name: google-adk-2026
description: Build code-first agents with Google Agent Development Kit (ADK). Use when you need explicit orchestration, evaluation, and Google ecosystem alignment without giving up model flexibility.
measurable_outcome: Select an ADK agent pattern, evaluation path, and deployment direction for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Google ADK (2026)

## Workflow

1. Read `references/sources.md` and confirm the current ADK quickstart, SDK surface, and A2A-related guidance.
2. Choose the right agent shape: single agent, multi-agent composition, or service exposed for interoperability.
3. Define tools, callbacks, and evaluation checkpoints before adding broad autonomy.
4. Decide whether the workload stays local, moves to Google-hosted infrastructure, or needs cross-framework interop through A2A.
5. Add a smoke test that exercises one tool call, one evaluation path, and one multi-agent handoff if applicable.

## When to Use

- You want code-first agent control with explicit orchestration.
- Google infrastructure or Gemini alignment matters.
- Multi-agent interoperability is on the roadmap.
- Evaluation should be part of the initial design, not an afterthought.

## Output Requirements

- Return the chosen ADK runtime pattern.
- State whether A2A is required now or later.
- Include one evaluation checkpoint.
- Include one deployment note covering local versus hosted execution.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
