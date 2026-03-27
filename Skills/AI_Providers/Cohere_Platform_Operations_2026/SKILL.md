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
name: cohere-platform-operations-2026
description: Integrate and operate Cohere APIs with current model and SDK guidance. Use when selecting Cohere models, implementing generation or embed/rerank pipelines, or planning migration across Cohere API updates.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Cohere Platform Operations (2026)

## Workflow

1. Check the current Cohere docs and SDK references in `references/sources.md`.
2. Choose the API surface that matches the workload: chat, embed, rerank, classify, or fine-tuning.
3. Prefer first-party SDKs and current request patterns before custom wrappers.
4. Define retries, quota assumptions, and fallback behavior for production traffic.
5. Validate with one representative request for the chosen pipeline before rollout.

## Output Requirements

- State the selected Cohere model family or API path.
- State the SDK choice and why.
- Include one operational guardrail tied to rate limits, retries, or migration risk.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
