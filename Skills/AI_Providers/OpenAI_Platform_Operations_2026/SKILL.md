---
name: openai-platform-operations-2026
description: Operate and integrate OpenAI APIs with current model, modality, and migration workflows. Use when building or updating OpenAI-based chat, reasoning, transcription, speech, image, or video pipelines, and when replacing deprecated models.
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
