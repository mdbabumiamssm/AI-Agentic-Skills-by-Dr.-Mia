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
name: anthropic-claude-operations-2026
description: Integrate and operate Anthropic Claude APIs with current model lifecycle guidance. Use when implementing Claude-based assistants, tool use, long-context reasoning, or when planning model upgrades based on release notes and deprecations.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Anthropic Claude Operations (2026)

## Workflow

1. Check release notes and deprecation notes in `references/sources.md`.
2. Choose model tier based on latency, cost, and reasoning depth.
3. Align tool-use strategy with official patterns before custom wrappers.
4. Define retries, safety checks, and fallback behavior.
5. Confirm expected behavior with a focused prompt/tool test.

## Output Requirements

- State chosen Claude model and why.
- State any known deprecation risk and deadline.
- List one fallback path for service continuity.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->