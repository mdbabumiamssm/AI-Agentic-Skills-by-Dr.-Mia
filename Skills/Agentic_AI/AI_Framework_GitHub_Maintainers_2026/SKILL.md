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
name: ai-framework-github-maintainers-2026
description: Track maintenance health of major AI framework repositories and choose reliable dependencies. Use when selecting orchestration stacks, planning upgrades, or reducing integration risk across fast-moving OSS ecosystems.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# AI Framework GitHub Maintainers (2026)

## Workflow

1. Start from the official repository list in `references/frameworks.md`.
2. Measure release recency, maintainer throughput, issue triage, and dependency churn before recommending a framework.
3. Prefer provider-backed or highly active first-party integrations over thin wrappers.
4. Capture framework-to-provider compatibility notes and known migration risks before upgrades.
5. Define smoke tests and rollback criteria for the core workflows that depend on the framework.

## Output Requirements

- List selected frameworks with a one-line maintenance rationale.
- Flag stale or risky dependencies and provide safer alternatives.
- Include a minimal upgrade test matrix for the main agent flows.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
