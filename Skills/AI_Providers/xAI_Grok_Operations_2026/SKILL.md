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
name: xai-grok-operations-2026
description: Integrate and operate xAI Grok APIs with current documentation and SDK guidance. Use when implementing xAI model calls, selecting Grok variants, or planning safe upgrades around API changes.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# xAI Grok Operations (2026)

## Workflow

1. Confirm the active docs and example repositories in `references/sources.md`.
2. Choose the Grok path that matches the workload: chat, vision, or voice.
3. Prefer official request examples and current endpoint patterns before custom wrappers.
4. Define retries, timeout budgets, and moderation or fallback expectations.
5. Run one representative smoke test on the target endpoint before broader rollout.

## Output Requirements

- State the selected Grok model family or endpoint path.
- State the official integration path being followed.
- Include one operational guardrail for rollout, retries, or model fallback.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
