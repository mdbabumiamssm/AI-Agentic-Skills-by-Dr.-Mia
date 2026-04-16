---
name: deepseek-reasoning-transition-agent
description: Reference older DeepSeek R1-specific reasoning examples and migrate them to current DeepSeek API operations guidance. Use when preserving an open-model reasoning example without treating raw chain-of-thought exposure as a stable product contract.
measurable_outcome: Reframe a DeepSeek R1 example as a current provider-aligned reasoning workflow with one safety or compatibility note.
allowed-tools:
  - read_file
  - run_shell_command
---

# DeepSeek Reasoning Transition Agent

This directory is a transition reference for model-specific reasoning examples.

## Workflow

1. Use `Skills/AI_Providers/DeepSeek_API_Operations_2026/` for current API, compatibility, and rollout guidance.
2. Preserve the useful part of the example: strong reasoning and code-oriented problem solving.
3. Do not assume that exposing raw reasoning traces is a durable or desirable application contract.
4. Validate structured outputs, retries, and compatibility against the current DeepSeek API docs before production use.

## Output Requirements

- State the provider-current DeepSeek skill that should anchor the implementation.
- Name one reasoning workload where the example still helps.
- Include one caution about exposing raw reasoning traces or relying on model-specific behavior.
