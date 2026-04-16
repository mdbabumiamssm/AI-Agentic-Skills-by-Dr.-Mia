---
name: gemini-flash-multimodal-transition-agent
description: Reference older Gemini Flash multimodal examples and migrate them to current Gemini operations guidance. Use when a pinned Gemini 2.0 example needs to be translated into the active Google GenAI SDK and model-selection path.
measurable_outcome: Translate a Gemini Flash example into a current multimodal implementation plan with one compatibility note.
allowed-tools:
  - read_file
  - run_shell_command
---

# Gemini Flash Multimodal Transition Agent

This directory is a model-specific reference example, not the main operational entry point.

## Workflow

1. Use `Skills/AI_Providers/Google_Gemini_Operations_2026/` for current Gemini model and SDK guidance.
2. Preserve the useful multimodal pattern: low-latency text, image, audio, or video handling.
3. Replace version-pinned assumptions with current model and SDK choices before shipping.
4. Run one text path and one multimodal path after migration.

## Output Requirements

- State which provider-current Gemini skill should anchor the implementation.
- Keep one useful multimodal pattern from the example.
- Include one compatibility note for moving from pinned examples to current SDK usage.
