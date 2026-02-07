---
name: google-gemini-operations-2026
description: Build and maintain Google Gemini API workflows using current docs and model catalog. Use when selecting Gemini models, implementing multimodal calls, or migrating code to the latest Google GenAI SDK patterns.
---

# Google Gemini Operations (2026)

## Workflow

1. Confirm active models and limits from official docs in `references/sources.md`.
2. Select Gemini model by context window, modality, and latency constraints.
3. Use official SDK initialization and request patterns.
4. Add schema validation for structured outputs.
5. Run a smoke test for text and one multimodal path when applicable.

## Output Requirements

- Provide chosen Gemini model and SDK path.
- Provide one compatibility note for previous integrations.
- Provide one operational guardrail (timeouts, retries, or quotas).
