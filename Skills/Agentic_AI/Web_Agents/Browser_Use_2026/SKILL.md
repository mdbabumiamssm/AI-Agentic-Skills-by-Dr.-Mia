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
name: browser-use-2026
description: Build browser-native AI task automation with Browser Use. Use when real web execution, persistent browser workspaces, and managed browser infrastructure are core to the workflow rather than incidental.
measurable_outcome: Select a Browser Use execution pattern, security boundary, and validation path for a concrete browser-agent workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Browser Use (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Browser Use SDK, workspace, and browser-agent execution model.
2. Decide whether the workload needs the high-level agent runtime or raw browser control through the browser surface.
3. Record the browser security boundary early: allowed domains, credentials, workspace files, proxies, and retention assumptions.
4. Treat browser execution as a substrate, not a replacement for broader orchestration when other frameworks already manage planning or approvals.
5. Add a smoke test that covers one end-to-end task, one file or workspace interaction, and one failure or auth boundary.

## When to Use

- Live browser execution is the core requirement.
- The workload needs persistent browser workspaces or managed infrastructure.
- Browser automation must be paired with explicit security and approval assumptions.
- You need a browser-native execution layer under a broader agent system.

## Output Requirements

- Return whether the design uses the agent surface or browser surface.
- State one auth, domain, or workspace boundary.
- Include one validation step for browser execution.
- Include one benchmark or reliability note for rollout context.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
