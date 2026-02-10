---
name: radgpt-radiology-reporter
description: Radiology Reporter
keywords:
  - radiology
  - report-generation
  - patient-friendly
  - summarization
  - explanation
measurable_outcome: Generate a patient-friendly explanation of a radiology report with <1% hallucination rate within 30 seconds.
license: MIT
metadata:
  author: Stanford Medicine
  version: "1.0.0"
compatibility:
  - system: Python 3.9+
allowed-tools:
  - run_shell_command
  - read_file
---

# RadGPT (Radiology Report Assistant)

An LLM-based agent designed to summarize and explain complex radiology reports for patients and clinicians.

## When to Use

*   **Patient Communication**: Converting technical findings into plain language.
*   **Clinician Review**: Highlighting critical findings (e.g., "Pneumothorax detected").
*   **Follow-up**: Suggesting appropriate next steps based on findings.

## Core Capabilities

1.  **Simplification**: Translates "bilateral opacity" to "cloudiness in both lungs".
2.  **Entity Extraction**: Identifies key anatomical structures and pathologies.
3.  **Q&A**: Answers follow-up questions about the report.

## Workflow

1.  **Input**: Raw text of the radiology report.
2.  **Process**: LLM summarizes and identifies key findings.
3.  **Output**: Structured summary or conversational explanation.

## Example Usage

**User**: "Explain this chest X-ray report to the patient."

**Agent Action**:
```bash
python -m radgpt.explain --report ./report.txt --target_audience patient
```
