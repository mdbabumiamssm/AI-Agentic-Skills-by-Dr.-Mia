---
name: anthropic-claude-operations-2026
description: Integrate and operate Anthropic Claude APIs with current model lifecycle guidance. Use when implementing Claude-based assistants, tool use, long-context reasoning, or when planning model upgrades based on release notes and deprecations.
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
