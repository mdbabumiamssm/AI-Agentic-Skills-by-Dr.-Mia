# Tutorial: Anthropic Co-Worker Workflows

**Stack:** Anthropic (Claude 3.5 Sonnet/Opus)
**Focus:** Complex Reasoning & Auditability
**Key Feature:** XML Tags & Thinking Blocks

## Introduction

In regulatory affairs, **getting the right answer isn't enough**. You need to know *why* the AI decided what it did. A black box cannot submit a drug application.

Anthropic's models are trained to be "Co-workers". They excel at:
1.  **Thinking Out Loud:** generating a `<thinking>` block before the final answer.
2.  **XML Protocol:** Using strict XML tags (e.g., `<citation>`) to structure documents.

## The Pattern: Regulatory Drafter

In `regulatory_drafter.py`, we simulate a "Co-worker" workflow:

### 1. The "Thinking" Phase
Before writing the draft, the agent outputs:
```xml
<thinking>
1. Analyze the Regulation (21 CFR).
2. Check the demographics in the clinical data.
3. Identify the discrepancy (Disease is Adult-only).
4. Decide to request a Waiver.
</thinking>
```
**Why this matters:** A human reviewer can read this. If the logic is flawed (e.g., "Disease is rare" vs "Disease is adult-only"), the human catches it *before* reading the draft.

### 2. The "Drafting" Phase
The agent produces the final text, using formal language and citations.

## How to use in Production

When prompting Claude, explicitly ask for this structure:

> "You are a Regulatory Affairs specialist. First, think through the strategy in `<thinking>` tags. Then, draft the submission in `<draft>` tags."

## Next Steps
1.  Get an Anthropic API Key.
2.  Install the SDK: `pip install anthropic`.
3.  Modify `regulatory_drafter.py` to accept a PDF text (the full FDA guidance) and see if it can extract specific constraints.
