<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

---
name: google-gemini-operations-2026
description: Build and maintain Google Gemini API workflows using current docs and model catalog. Use when selecting Gemini models, implementing multimodal calls, or migrating code to the latest Google GenAI SDK patterns.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Google Gemini Operations (2026)

## Workflow

1. Confirm active models and limits from official docs in `references/sources.md`.
2. Select Gemini model by context window, modality, and latency constraints.
3. Use official SDK initialization and request patterns.
4. Add schema validation for structured outputs.
5. Run a smoke test for text and one multimodal path when applicable.

## Output Requirements

- Provide chosen Gemini model and SDK path.
- Provide one compatibility note for previous integrations.
- Provide one operational guardrail (timeouts, retries, or quotas).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->