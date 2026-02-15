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
name: ai-provider-github-maintainers-2026
description: Track official AI provider GitHub repositories and release health for reliable integrations. Use when auditing SDK freshness, selecting maintained repos, or planning upgrades based on release cadence and issue velocity.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# AI Provider GitHub Maintainers (2026)

## Workflow

1. Start from the official repository list in `references/providers.md`.
2. Check release cadence and open issue trends for each dependency.
3. Prefer provider-maintained SDKs over third-party wrappers.
4. Pin versions for production, then test upgrade candidates in staging.
5. Record upgrade windows and rollback procedures.

## Output Requirements

- List selected repositories with short rationale.
- Flag stale repositories and replacement options.
- Include an upgrade test checklist.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->