---
name: claude-reasoning-transition-agent
description: Reference older Anthropic model-specific reasoning patterns and migrate them to the current Anthropic operations skill. Use when maintaining a legacy Claude 3.7-specific example or translating it to current Claude model selection guidance.
measurable_outcome: Convert a model-pinned Claude reasoning example into a provider-current implementation plan with one migration note.
allowed-tools:
  - read_file
  - run_shell_command
---

# Claude Reasoning Transition Agent

This directory is retained as a model-specific reference example. It should not be the primary starting point for current Anthropic integration work.

## Workflow

1. Use `Skills/AI_Providers/Anthropic_Claude_Operations_2026/` for active model and deprecation guidance.
2. Keep the useful pattern, not the pinned model assumption: long-form reasoning, coding assistance, and tool-aware analysis.
3. Remove any hard dependency on a single Claude model line before production rollout.
4. Validate the migrated workflow against the current Anthropic model overview and release notes.

## Output Requirements

- State the current Anthropic skill that should replace this model-pinned reference.
- Keep one useful reasoning pattern from the example.
- Include one migration note about model pinning, tool use, or latency.
