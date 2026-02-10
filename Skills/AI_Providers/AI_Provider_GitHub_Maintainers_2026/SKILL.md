---
name: ai-provider-github-maintainers-2026
description: Track official AI provider GitHub repositories and release health for reliable integrations. Use when auditing SDK freshness, selecting maintained repos, or planning upgrades based on release cadence and issue velocity.
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
