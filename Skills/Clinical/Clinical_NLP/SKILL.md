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

---
name: 'clinical-nlp-extractor'
description: 'Extracts medical entities (Diseases, Medications, Procedures) and patient outcome mentions with temporal anchoring from unstructured clinical text using regex and simple rules (or LLM wrappers).'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---


# Clinical NLP Entity Extractor

The **Clinical NLP Skill** converts free-text clinical notes into structured data. It identifies key medical entities like problems/diagnoses, medications, and procedures.

## When to Use This Skill

*   When analyzing unstructured EHR notes.
*   To populate a patient's problem list or medication reconciliation.
*   To de-identify text (phi-removal) - *Basic version*.

## Core Capabilities

1.  **NER (Named Entity Recognition)**: Extracts Problems, Drugs, Procedures.
2.  **Negation Detection**: (Basic) Checks if a finding is denied ("No fever").
3.  **Structuring**: Returns JSON format compatible with FHIR/USDL.
4.  **Incidentaloma Identification**: Detect incidental imaging findings requiring follow-up across multiple anatomies in radiology reports. Park et al. (J Biomed Inform, Apr 2026, PMID 42061667) compare LLM-based vs supervised approaches — useful exemplar when scoping radiology-report NLP pipelines.
5.  **Outcome Identification with Temporal Anchoring**: Detect outcome mentions in clinician notes, normalize timing, distinguish documented events from speculation, and return evidence spans for audit.
6.  **Clinician-Note Outcome Timing Extraction**: Identify patient outcome mentions in clinician notes, normalize event timing, distinguish historical/current/future events, link extracted spans to structured variables, and validate against adjudicated chart-review labels.
7.  **Registry-Oriented Rare-Cancer RWD Extraction**: For low-prevalence cancer cohorts such as bone sarcoma, define the registry schema and field-level data dictionary before extracting and normalizing diagnoses, pathology type, tumor size, localization, grade, primary resection, treatment timelines, and outcomes; link each note-to-field value to its evidence span and source-note identifier; normalize documented dates and relative temporal expressions; distinguish absent documentation, unknown values, and explicit negative findings; route ambiguous or sampled values for human validation; monitor extraction quality with field-level and aggregate audit metrics such as precision and recall without assuming preset benchmarks; and preserve temporal anchors, missingness states, and validation status as provenance.
8.  **Oncology Consultation-Document Survival Modeling**: For prognostic modeling from initial oncology consultation documents, compare zero-shot prompting with task-specific fine-tuning before selecting an approach; use censoring-aware survival labels, leakage checks for post-baseline or outcome-revealing text, held-out calibration assessment, and clear limits that outputs are investigational decision-support signals rather than direct clinical-use recommendations.
9.  **Longitudinal Outcome-Event Identification and Timing**: Extract explicit and inferred outcome events with negation status; anchor onset and resolution as dates or intervals when timing is uncertain; separate note time from event time; reconcile, update, and deduplicate events across longitudinal notes; retain clinician-auditable evidence spans; and evaluate event-level precision and recall without assuming preset benchmarks.
10. **Rare-Cancer Longitudinal Registry Abstraction**: Extract registry fields from longitudinal rare-cancer notes using ontology-backed field schemas; resolve event chronology across note and event dates; encode explicit negatives, unknown values, and absent documentation as distinct missingness states; assign field-level confidence with supporting evidence spans; route uncertain records for clinician validation; and compare extraction results against manual abstraction without assuming preset performance benchmarks.
11. **Bone-Sarcoma Registry Extraction Pattern**: For rare-cancer registries, use bone sarcoma as a worked pattern: define the abstraction schema before extraction, link temporal events across records, normalize pathology and treatment fields, attach uncertainty flags, validate behavior for low-prevalence data, and manually audit a sample of extracted records.
12. **LLM Rare-Neoplasm RWD Extraction Workflow**: Use LLMs to structure oncology notes for rare-neoplasm real-world-data collection; for bone sarcoma, capture registry-defined disease-specific variables with source-note identifiers and evidence spans for each extracted field; preserve prompt, model, run, and reviewer provenance where available; validate extracted values against manual abstraction before registry loading; and quantify field-level extraction error so uncertain or high-error variables can be corrected or withheld from registry use.

## Workflow

1.  **Input**: A string of clinical text or a text file.
2.  **Process**: Tokenizes and matches against patterns/dictionaries.
3.  **Output**: JSON list of entities with spans and types.

## Example Usage

**User**: "Extract entities from this note."

**Agent Action**:
```bash
python3 Skills/Clinical/Clinical_NLP/entity_extractor.py \
    --text "Patient has diabetes type 2. Prescribed Metformin 500mg. No chest pain." \
    --output entities.json
```

## References

- Park N, Ahmed F, Sun Z, Lybarger K, Breinhorst E. Automated identification of incidentalomas requiring follow-up: A multi-anatomy evaluation of LLM-based and supervised approaches. J Biomed Inform, 2026 Apr 28. https://pubmed.ncbi.nlm.nih.gov/42061667/
- Abdullahi T, Hamzeh A, Sears I, Abadi N, Singh R. Identifying and timing patient outcomes in clinician notes using large language models. Artif Intell Med, 2026 Jun. https://pubmed.ncbi.nlm.nih.gov/41886942/
- Teterycz P, Rynkun S, Szostakowski B, Wągrodzki M, Rutkowski P. Accelerating real-world data collection using large language models in rare neoplasms: a bone sarcoma example. ESMO Real World Data Digit Oncol, 2026 Jun. https://pubmed.ncbi.nlm.nih.gov/42021926/
- Phaterpekar T, Zeng Z, Mali Y, Leung B, Ho C. Investigating fine-tuning versus zero-shot learning for general large language models when predicting cancer survival from initial oncology consultation documents. ESMO Real World Data Digit Oncol, 2026 Jun. https://pubmed.ncbi.nlm.nih.gov/42004490/

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
