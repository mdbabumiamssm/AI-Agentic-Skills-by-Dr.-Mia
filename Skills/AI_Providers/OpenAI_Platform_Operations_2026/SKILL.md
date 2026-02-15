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
name: openai-platform-operations-2026
description: Operate and integrate OpenAI APIs with current model, modality, and migration workflows. Use when building or updating OpenAI-based chat, reasoning, transcription, speech, image, or video pipelines, and when replacing deprecated models.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# OpenAI Platform Operations (2026)

## Workflow

1. Read `references/sources.md` and verify the current model and endpoint status before coding.
2. Select API path by workload: text/reasoning, audio, image, or video.
3. Prefer official SDK examples, then adapt for local constraints.
4. Add explicit migration notes when replacing old or deprecated models.
5. Validate with a minimal end-to-end request before broader refactors.

## Output Requirements

- Return selected model and endpoint choice with a one-line rationale.
- Include rate-limit, timeout, and retry assumptions.
- Include a short migration note if any model replacement was made.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->