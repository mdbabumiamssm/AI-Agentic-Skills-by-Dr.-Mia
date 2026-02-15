<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Tutorial: Building OpenAI Health Agents with Structured Outputs

**Stack:** OpenAI (GPT-4o)
**Focus:** Reliability & Compliance
**Key Feature:** JSON Schema Enforcement

## Introduction

In healthcare, "hallucinations" are dangerous. You cannot have an AI suggest a made-up drug or vague advice.

OpenAI's **Structured Outputs** (JSON Mode) solve this by forcing the model to adhere to a strict blueprint. If the model tries to output text that doesn't fit your schema, the API request fails (or is corrected).

## The Pattern: Care Copilot

In `care_copilot.py`, we demonstrate the standard pattern for clinical intake:

1.  **Define the Schema:** We create a precise definition of what a "Triage Result" looks like.
2.  **Strict Mode:** We set `"strict": True`. This guarantees the output matches the schema 100%.
3.  **Process:** We ingest messy natural language ("my tummy hurts") and get clean, database-ready JSON.

## Code Walkthrough

```python
triage_schema = {
    "type": "object",
    "properties": {
        "priority": {"enum": ["LOW", "HIGH"]}, # Constrained Choice
        "action": {"type": "string"}
    },
    "required": ["priority", "action"],
    "additionalProperties": False # No hallucinations allowed
}
```

## Why this is "Enterprise Grade"

*   **Type Safety:** The output plugs directly into your EMR/EHR system.
*   **Guardrails:** You can restrict the "Priority" field to specific enums, preventing the AI from inventing new triage levels.
*   **Auditability:** Every decision is structured and logged.

## Next Steps

1.  Get an OpenAI API Key.
2.  Uncomment the `client = OpenAI(...)` lines in `care_copilot.py`.
3.  Expand the schema to include `vitals` (heart rate, bp) extraction.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->