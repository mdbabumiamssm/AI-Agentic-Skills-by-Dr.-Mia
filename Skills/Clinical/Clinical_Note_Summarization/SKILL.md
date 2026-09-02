---
name: clinical-note-summarization
description: Structure raw clinical notes into SOAP-format summaries with explicit contradictions, missing data, and ICD-linked assessments using the provided prompt + usage script.
measurable_outcome: Produce SOAP markdown and JSON outputs covering all four sections with at least 95% note coverage and explicit missing information within 2 minutes per note.
allowed-tools:
  - read_file
  - run_shell_command
---

<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

## At-a-Glance
- **description (10-20 chars):** SOAP builder
- **keywords:** clinical-notes, SOAP, guardrails, ICD10, gaps
- **measurable_outcome:** Produce SOAP markdown + JSON (when requested) covering all four sections with ≥95% note coverage and explicit missing info in ≤2 minutes per note.

## Inputs
- `note_text` (dictation, OCR, or EHR export) and optional `patient_context` metadata.
- `output_format` (`markdown` default, `json` when downstream validators need schema).

## Outputs
1. Structured SOAP summary with Subjective/Objective/Assessment/Plan bulleting.
2. Alerts plus missing-information checklist.
3. Optional JSON payload using schema from README.

## Core Capabilities
- Evaluate AI-generated clinical note quality for factual consistency, omitted clinically relevant facts, temporal accuracy, diagnosis/medication/procedure fidelity, note completeness, and rubric-based human review; treat automated benchmarks as supportive, caveated quality checks rather than definitive clinical validation.
- For AI-assisted urologic documentation, maintain source-to-note traceability; explicitly check for omissions and contradictions; detect and flag templated or copied-forward text; require responsible clinician review, editing, and attestation; and limit generated content to documentation support rather than autonomous clinical judgment.
- For generated or transformed urologic documentation, preserve procedure- and operative-note structure; verify staging and laterality; reconcile medications and devices; detect and flag contradictions; attribute statements to their source documentation; and require responsible clinician review and sign-off before record entry or patient-facing use.
- For urologic workflows, preserve procedural details, laterality, staging, and device information; keep drafts linked to source documentation; detect and flag omissions and copy-forward errors; and require clinician sign-off before medico-legal content is finalized.
- Apply benchmark-driven quality evaluation for AI-generated notes across correctness, omissions, factual consistency, note completeness, risk-of-harm scoring, and human review rubric selection.
- Compare automated note-quality checks against clinician review in benchmark-style evaluation across factual consistency, omission detection, note completeness, safety/harm review, and rubric-based scoring.
- For urologic consultation, procedure, and follow-up notes, check factual completeness against the source, detect and flag unsupported content, require clinician review and sign-off before record entry, and prohibit autonomous documentation without human oversight.
- Support human-reviewed AI drafting and summarization for urologic documentation, including consultation, procedure, and follow-up notes, while requiring verification before record entry.
- Structure LLM-assisted urologic note drafts with explicit history and procedure fields, preserve specialty vocabulary, enforce privacy constraints for source text and generated drafts, require clinician sign-off, and evaluate draft quality against comparable human-authored documentation.
- For urologic specialty documentation, preserve procedure-specific terminology; distinguish summarization from autonomous documentation; validate omissions and unsupported additions; and require clinician sign-off before record entry.
- For AI-assisted urologic operative notes, clinic documentation, and patient instructions, require human review for specialty-specific terminology, omitted negatives, billing/coding risk, and medicolegal traceability before record entry or patient release.
- For AI-drafted urologic clinical notes, procedure documentation, and coding-support language, require human review checkpoints that validate specialty terminology and source support, with explicit checks for omitted exam findings, copied-forward content, and overconfident generated plans before record entry or coding use.
- Apply specialty-documentation safety checks for urologic note drafts, including clinician oversight, source-note traceability, procedure-specific omission review, and separation from billing or medico-legal finalization.
- Use AI-assisted urologic documentation only as structured note drafting support that may reduce documentation burden: preserve visit-specific terminology across consultation, procedure, follow-up, operative, clinic, and patient-instruction contexts; complete privacy review for source text and generated drafts; and prohibit using generated notes without clinician verification, editing, and sign-off.
- For urologic documentation use cases, perform specialty terminology checks, enforce procedure-note constraints, draft with PHI-safe handling, maintain audit trails for source-to-output review and edits, and require responsible clinician sign-off before generated documentation enters the record.
- For LLM-assisted urologic documentation, align drafts to the applicable templated note sections; preserve procedure-specific terminology; check required sections for omissions; safeguard PHI in source text, prompts, and generated drafts; and require clinician oversight, editing, and sign-off before record entry.

## Workflow
1. **Load system prompt:** `prompt.md` enforces no hallucinations + data gap surfacing.
2. **Normalize input:** Pre-clean vitals, labs, and timeline context when available.
3. **Generate summary:** Call preferred LLM (OpenAI, Anthropic, Gemini, OSS) using `usage.py` as a template.
4. **Validate:** Cross-check extracted values vs. source text and ensure contradictions/missing data are spelled out.
5. **Deliver output:** Provide markdown + JSON as required and log PHI handling steps.

## Guardrails
- Never invent findings; state "not provided" explicitly.
- Mark outputs as documentation support only—not clinical decisions.
- Strip/re-mask PHI before storing prompts/responses.

## Human Oversight Before Record Entry
- **Factual verification:** Compare symptoms, history, examination findings, test results, diagnoses, medications, procedures, and dates against the source documentation.
- **Unsupported additions:** Remove or explicitly flag statements, inferred negatives, or recommendations that are not supported by the source.
- **Coding implications:** Require clinician or qualified coding review of diagnosis and procedure wording and any resulting coding implications.
- **Privacy:** Confirm that PHI handling, storage, transmission, and access follow applicable organizational privacy controls.
- **Clinician attestation:** Require the responsible clinician to review, edit as needed, and attest the final documentation.
- **Failure handling:** Withhold the generated note from the record when verification fails, source data are unavailable, or the system errors; route the case to manual documentation and correction.

## Specialty-Documentation Safety
- Require clinician oversight for urologic documentation drafts before record entry or patient-facing use.
- Preserve traceability to source notes for diagnoses, procedures, findings, and follow-up plans.
- Check omissions of procedure-specific details and specialty terminology before clinician review.
- Keep AI-generated drafts separate from billing or medico-legal finalization until clinician attestation is complete.

## References
- For detailed schema, guardrails, and integration snippets see `README.md`, `prompt.md`, and `usage.py`.
- PubMed PMID 41955894: https://pubmed.ncbi.nlm.nih.gov/41955894/
- PubMed PMID 42067659: https://pubmed.ncbi.nlm.nih.gov/42067659/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
