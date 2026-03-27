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
name: mistral-platform-operations-2026
description: Integrate and operate Mistral APIs with current model catalog and SDK guidance. Use when implementing Mistral chat or embedding workflows, selecting models, or migrating between API and SDK versions.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Mistral Platform Operations (2026)

## Workflow

1. Confirm the current model catalog and SDK guidance in `references/sources.md`.
2. Choose the Mistral capability set that matches the workload: chat, reasoning, embeddings, coding, or document AI.
3. Prefer first-party SDK patterns for structured outputs, citations, and tool use.
4. Define retry behavior, timeout budgets, and rollout notes before model changes.
5. Run one focused smoke test against the chosen model path before wider changes.

## Output Requirements

- State the selected Mistral model or API path.
- State the SDK path and any compatibility note.
- Include one operational guardrail for retries, quotas, or fallback.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
