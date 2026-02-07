---
name: frontier-oss-models-2026
description: Work with leading open or openly available model ecosystems from Meta, Mistral, DeepSeek, and xAI. Use when evaluating self-hosted or hybrid deployment paths, model tradeoffs, and integration patterns across non-single-vendor stacks.
---

# Frontier OSS Models (2026)

## Workflow

1. Identify target deployment mode: hosted API, self-hosted, or hybrid.
2. Review provider-specific docs and repos in `references/sources.md`.
3. Compare model choice by licensing, latency, and hardware profile.
4. Define eval set before changing production model baselines.
5. Roll out through canary traffic with measurable success thresholds.

## Output Requirements

- State selected model family and deployment mode.
- List one licensing/compliance check.
- List one rollback condition tied to quality or cost.
