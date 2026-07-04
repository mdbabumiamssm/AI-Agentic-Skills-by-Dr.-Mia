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



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'rare-neoplasm-rwd-llm-extraction'
description: 'Extract registry-ready rare-neoplasm variables from clinical text using schema-first LLM workflows with temporal normalization, ontology mapping, provenance, adjudication, privacy, and validation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Rare Neoplasm Real-World Data LLM Extraction

## Overview

Extract structured real-world data from oncology notes, pathology reports, operative records, imaging reports, and treatment summaries for rare neoplasms. Use a schema-first workflow to handle sparse documentation and terminology variation while preserving evidence provenance, uncertainty, privacy, and a clear route to human review.

## When to Use This Skill

- Build or update a rare-cancer registry from unstructured clinical documents.
- Extract bone sarcoma or other rare-neoplasm diagnoses, pathology, staging, treatment, response, recurrence, and outcomes.
- Harmonize institution-specific tumor terminology with oncology ontologies.
- Reconstruct longitudinal cancer timelines from fragmented records.
- Audit an LLM extraction pipeline against manually curated registry data.
- Produce review queues for ambiguous, conflicting, or low-confidence variables.
- Do not use the workflow as a substitute for clinical diagnosis, treatment decisions, or source-document verification.

## Core Capabilities

1. **Define the extraction contract.** Specify the cohort, tumor family, source types, index event, field definitions, cardinality, allowed values, units, and required evidence before processing records. Version the schema and ontology releases so outputs remain reproducible.

2. **Apply privacy controls.** Process only authorized records in an approved environment, enforce least-privilege access, minimize retained identifiers, avoid placing protected health information in prompts or logs beyond operational need, and follow applicable institutional and legal requirements. De-identify or pseudonymize exports unless identifiable data are explicitly required and authorized.

3. **Prepare source documents.** Preserve document identifiers, document type, authoring service, encounter date, specimen date, and original text offsets. Remove duplicate documents, identify amended reports, and retain report versions so later corrections can supersede earlier findings without erasing provenance.

4. **Extract schema-bound facts.** Return only requested fields in a machine-valid structure such as JSON or a registry table. For each assertion, capture the raw value, normalized value, status, source document, evidence span, character or line offsets when available, and extraction timestamp.

5. **Link every value to evidence.** Use the shortest sufficient verbatim span supporting each extracted fact. Distinguish direct statements from inferred values; do not present an inference as explicitly documented. If a derived field is required, record the derivation rule and all supporting source fields.

6. **Normalize time and clinical episodes.** Convert dates to ISO 8601 when the source supports it, retain the original expression, and represent partial dates without inventing missing precision. Anchor relative expressions to a documented reference date, associate events with the correct tumor episode, and flag impossible or conflicting sequences.

7. **Map rare-tumor terminology.** Preserve the original disease wording and map it to the project-approved ontology, such as ICD-O, NCI Thesaurus, or SNOMED CT, including code, label, ontology version, and mapping confidence. Use morphology, primary site, molecular findings, and pathology context; abstain rather than force an unsupported code.

8. **Resolve negation, uncertainty, and context.** Separate confirmed findings from suspected, ruled-out, historical, family-history, planned, canceled, and hypothetical statements. Track whether a statement applies to the patient, a specific specimen, a metastatic site, or another tumor episode.

9. **Represent missingness explicitly.** Distinguish `not_documented`, `not_assessed`, `not_applicable`, `unknown`, `redacted`, `conflicting`, and `extraction_failed`. Never convert absent text into a negative clinical finding or silently impute a value.

10. **Calibrate confidence and abstention.** Assign field-level confidence using documented evidence quality, terminology fit, temporal consistency, and cross-document agreement. Configure thresholds by field risk: auto-accept high-confidence facts, route intermediate-confidence facts to review, and abstain on unsupported or contradictory facts.

11. **Reconcile longitudinal conflicts.** Prefer final pathology over preliminary interpretations when appropriate, amended reports over superseded versions, and more specific evidence over generic summaries. Preserve all conflicting assertions and record the deterministic precedence rule or human decision used to select the registry value.

12. **Require human adjudication where needed.** Create a review queue containing the candidate value, alternatives, evidence spans, source metadata, conflict reason, and proposed resolution. Require qualified review for low-confidence mappings, discordant pathology, uncertain primary site, duplicate tumor episodes, and safety-critical fields.

13. **Validate against curated data.** Compare outputs with an independently manually curated registry using a locked schema and adjudication guide. Report field-level completeness, exact or clinically meaningful agreement, precision, recall, error categories, abstention rate, and review burden as applicable; stratify by document type and tumor subtype without inventing performance claims.

14. **Use the bone sarcoma rare-neoplasm example as a validation case.** Validate schema-first variable extraction for rare-neoplasm real-world data collection with a bone sarcoma scenario that exercises sarcoma-specific schemas, temporality and event ordering, ontology normalization, evidence provenance capture, registry readiness, adjudication workflow routing, and field-level metrics for rare-neoplasm chart abstraction, aligned to the 2026 PubMed-indexed example without inventing unsupported benchmarks or claims.

15. **Use bone sarcoma RWD extraction examples.** For sparse rare-neoplasm records, capture longitudinal outcome variables only when documented, map disease, treatment, and event terms to approved ontologies, normalize partial or relative dates, link each abstraction to source documents and evidence spans, and route ambiguous or conflicting registry-ready variables to human adjudication.

16. **Maintain an audit trail.** Record model and prompt versions, schema and ontology versions, preprocessing steps, source hashes or stable identifiers, extraction status, reviewer actions, and final disposition. Make each final value traceable to source evidence and any transformation or adjudication.

17. **Run release checks.** Reject outputs that fail schema validation, contain unsupported values, lack required provenance, expose unauthorized identifiers, or violate chronology constraints. Sample accepted, rejected, and abstained cases for periodic human quality review and monitor drift when data sources or models change.

## Inputs / Outputs

### Inputs

- Authorized clinical text or document references, with stable document and patient or case identifiers.
- A versioned rare-neoplasm extraction schema and data dictionary.
- Cohort definition, index-date rules, tumor-episode rules, and document precedence rules.
- Approved ontology names, versions, code sets, and local terminology mappings.
- Field-specific confidence thresholds and human-review criteria.
- Privacy, retention, access-control, and export requirements.
- Optional manually curated registry records for validation.

### Outputs

- Schema-valid patient-, tumor-, specimen-, and event-level records.
- For each extracted field: raw text value, normalized value, evidence span, source identifier, temporal anchor, ontology mapping, confidence, and assertion status.
- Explicit missingness, abstention, conflict, and extraction-failure states.
- A chronological tumor and treatment event timeline with unresolved inconsistencies flagged.
- A human-adjudication queue with evidence and resolution history.
- Validation and error-analysis summaries against the curated reference set.
- A reproducibility and audit manifest containing model, prompt, schema, ontology, preprocessing, and review versions.

## References

- Teterycz P, Rynkun S, Szostakowski B, Wągrodzki M, Rutkowski P. “Accelerating real-world data collection using large language models in rare neoplasms: a bone sarcoma example.” *ESMO Real World Data and Digital Oncology*. 2026 Jun. PubMed: https://pubmed.ncbi.nlm.nih.gov/42021926/
