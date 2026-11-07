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
name: 'longitudinal-clinical-memory-evaluation'
description: 'Construct and evaluate patient-state memory across medical dialogue sessions, testing dated state changes, person attribution, and retrieval of safety-relevant evidence.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Longitudinal Clinical Memory Evaluation

## Overview

Build an auditable research representation of patient state across dialogue sessions and evaluate whether retrieval preserves time, person identity, and safety-relevant evidence. The supplied MedGraph finding motivates these three failure categories and graph-based state histories; the workflow below is an operational research adaptation, not a claim of exact reproduction or validated clinical deployment guidance.

## When to Use This Skill

- Comparing memory architectures for assistants that follow patients across visits.
- Auditing answers that reuse discontinued treatments or superseded clinical states.
- Investigating confusion between the patient, relatives, caregivers, or other people.
- Testing whether recorded safety facts reach the answering model consistently.
- Designing chronological retrieval tests, graph-memory ablations, or small reproducible pilots.

## Core Capabilities

1. **Define a bounded evaluation.** Select a small, explicit cohort or synthetic fixture set, session cutoff, questions, and comparison configuration. Use authorized, de-identified dialogue in an approved environment; do not send identifiable records to external tools. Reserve approximately 2 minutes for scope, 4 for extraction, 3 for retrieval, 4 for scoring, and 2 for reporting. If resources exceed this budget, report the completed subset and pending work without claiming full evaluation.

2. **Extract dated, attributable evidence.** Create one immutable record per asserted fact, preserving its original dialogue span. Distinguish the speaker from the person described; a caregiver speaking about the patient does not make the fact a caregiver fact. Assign unresolved subjects an explicit unknown value instead of defaulting to the patient. Preserve negation, uncertainty, hypothetical statements, and reported versus confirmed status.

   Minimum fact fields:

   ```text
   fact_id, patient_partition_id, subject_id, subject_relation, speaker_id,
   session_id, source_span_id, source_text, recorded_at,
   effective_start, effective_end, date_precision,
   concept_type, concept, value, unit, assertion_status,
   clinical_status, safety_tag, extraction_uncertainty
   ```

   Keep unknown dates null. Resolve relative dates only when the session provides a reliable anchor, and record the derivation. Separate when a statement was recorded from when its clinical event occurred. Do not silently normalize ambiguous doses, units, or names.

3. **Maintain state histories.** Partition records by patient and subject before linking medication, diagnosis, and measurement events. Connect records only when concept and identity matching are supported, using typed links such as `updates`, `discontinues`, `resolves`, `corrects`, or `contradicts`. Retain earlier records and their provenance. A later mention alone does not supersede an earlier fact; absence of a mention does not establish resolution. Measurements remain dated observations rather than indefinitely current values. Preserve conflicting assertions for review instead of inventing a reconciled state.

4. **Retrieve the appropriate state and safety context.** Represent each question with both a knowledge cutoff and the clinical time being asked about. Admit only facts recorded by the cutoff, then resolve state at the requested clinical time; later disclosures must not leak into earlier-session tests. Rank eligible patient chains as units, select supported states, and attach their source spans. Keep family-history evidence explicitly labeled with its actual subject.

   Assemble a separate profile of all recorded patient safety facts available at the cutoff, including their dates, status, uncertainty, and provenance. Examples of evaluation labels include allergies, adverse reactions, and explicitly recorded contraindications; these labels are not an exhaustive clinical checklist. Retain historical or disputed entries with those labels rather than presenting them as active. Include this profile with each question and audit the final model input for truncation. Missing documentation means unknown, not absence of risk. If the question time or evidence is insufficient, return a qualified answer or abstention.

5. **Create independent reference judgments.** Annotate expected subject, time-specific state, required supporting spans, required safety facts, and acceptable answer or abstention from the source dialogue before inspecting system outputs. Include medication changes, resolved or disputed diagnoses, repeat measurements, relatives with similar conditions, indirect speakers, safety facts from early sessions, uncertain dates, and late corrections. Label synthetic fixtures as synthetic. Keep evaluation questions and answer labels out of memory construction, and separate patients across tuning and held-out evaluation. Record annotator disagreement and adjudication status.

6. **Measure retrieval and answer failures separately.** For each question, save retrieved fact IDs, final context IDs, answer, supporting citations, reference labels, and error flags. Score the following operational metrics; they are proposed evaluation definitions, not claimed MedMemoryBench definitions.

   | Metric | Numerator / denominator |
   | --- | --- |
   | Stale-state retrieval rate | Time-sensitive questions whose final context presents a superseded state as applicable / all time-sensitive questions |
   | Wrong-person retrieval rate | Questions whose final context misattributes another person's fact to the patient / all evaluated questions |
   | Safety-evidence recall | Required safety fact instances present with correct subject and status in final context / all required safety fact instances |
   | Supporting-evidence recall | Required reference evidence instances present in final context / all required reference evidence instances |
   | Answer correctness | Answers meeting the predeclared subject, time, content, and evidence rubric / all scored questions |
   | Abstention rate | Abstained questions / all evaluated questions |

   Report counts and denominators, with `N/A` for zero denominators. A historical fact clearly labeled historical is not a stale-state error. Properly attributed family history is not a wrong-person error. Allow multiple failure flags per question, and additionally mark whether each retrieval failure propagated into the answer. Score abstentions against answerability labels; do not treat every abstention as correct. Distinguish extraction omissions, incorrect links, retrieval omissions, context truncation, and answer-generation errors.

7. **Compare and report reproducibly.** When resources permit, compare the graph approach with a declared baseline such as chronological dialogue retrieval, using the same available sessions, questions, answer model, and context budget. Ablate temporal resolution, subject filtering, and safety-profile inclusion separately to test their contributions. Record model/version, prompts, retrieval settings, token limits, seeds where supported, construction time, answering latency, and token usage. Report observed results only; mark unexecuted comparisons as pending. Small pilots support debugging, not general clinical effectiveness claims.

## Inputs / Outputs

**Inputs**

- Authorized dialogue sessions with patient partition keys, session timestamps, speaker identifiers, and stable source-span references.
- Questions with knowledge cutoffs and requested clinical times; explicit handling rules for ambiguous time references.
- Reference annotations or sufficient source evidence to create them, including person attribution and safety dependencies.
- Available memory/retrieval implementations, model configuration, context budget, and evaluation scope. If no model is available, perform a labeled manual retrieval audit and mark answer-generation metrics unmeasured.

**Outputs**

- A fact ledger and state-chain representation retaining timestamps, attribution, contradictions, and source provenance.
- Per-question safety profiles and retrieval traces showing exactly what entered the final answer context.
- A scored question table with expected and observed states, evidence coverage, abstentions, and independent failure flags.
- A concise run report with scope, configuration, metric counts, representative errors, completed comparisons, limitations, and next experiments.

Return these as structured tables or JSON in the response unless artifact files are explicitly requested. A valid pilot output contains traceable fact records, at least one scored test for each of the three target failure categories, explicit denominators, and a clear distinction between measured, manually assessed, and unexecuted results. If required inputs are missing, return an input-gap report and unexecuted evaluation plan rather than fabricated scores. Completion within 15 minutes describes the bounded workflow, not a guarantee of full benchmark execution or clinical readiness.

## References

- [MedGraph: Maintaining Patient State with Graph Memory for Long-Term Medical Dialogue](http://arxiv.org/abs/2611.08433v1), arXiv:2611.08433v1, submitted November 5, 2026. Source for the three memory failure categories and the proposed dated, subject-aware chain retrieval with a safety profile. Treat reported benchmark findings as research evidence; consult the full methods and available artifacts before claiming reproduction.
