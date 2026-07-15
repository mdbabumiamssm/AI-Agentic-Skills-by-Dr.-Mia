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
name: 'oncology-neurosymbolic-trial-matching'
description: 'Match oncology patients to clinical trials using knowledge-graph retrieval, symbolic eligibility reasoning, specialized agents, conflict resolution, and clinician-auditable evidence.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Oncology Neuro-Symbolic Trial Matching

## Overview

This skill produces explainable oncology clinical trial candidate matches by combining oncology-specific knowledge representation, symbolic eligibility evaluation, and coordinated specialist-agent review. It is grounded in a prospective evaluation involving 3,804 patients and is designed to support clinician review rather than replace clinical judgment.

## When to Use This Skill

- A clinician or trial-navigation team needs candidate trials for a patient with cancer.
- Eligibility criteria must be evaluated against structured and unstructured patient data.
- Trial retrieval must account for disease ontology, biomarkers, prior therapies, stage, and geography.
- Multiple agents or matching methods disagree and require explicit conflict resolution.
- A trial shortlist needs criterion-level evidence, uncertainty labels, and an audit trail.
- Prospective matching operations need a reproducible human-in-the-loop workflow.

## Core Capabilities

1. **Patient representation:** Normalize diagnosis, histology, stage, molecular findings, treatment history, performance status, organ function, comorbidities, age, and location while preserving source provenance and dates.
2. **Oncology knowledge-graph retrieval:** Connect patient concepts to cancer types, variants, therapies, trial interventions, sites, and eligibility concepts using ontology-aware expansion without treating inferred relationships as confirmed facts.
3. **Symbolic eligibility reasoning:** Translate inclusion and exclusion criteria into explicit rules and classify each criterion as met, not met, unknown, or not applicable.
4. **Specialist-agent decomposition:** Assign focused review roles for disease fit, biomarker fit, treatment history, laboratory and clinical constraints, temporal logic, and geographic feasibility.
5. **Evidence-grounded synthesis:** Require every match or exclusion claim to cite the relevant patient datum and trial criterion; never fill missing clinical facts by assumption.
6. **Conflict resolution:** Compare agent conclusions, expose contradictory interpretations, apply precedence to directly documented evidence, and route unresolved conflicts for clinician review.
7. **Ranking with uncertainty:** Rank candidates by eligibility confidence, clinical relevance, recruiting status, location, and unresolved-data burden without presenting the ranking as a treatment recommendation.
8. **Clinician review package:** Produce a compact, auditable report containing candidate trials, criterion-level decisions, missing-data requests, conflicts, provenance, and review timestamps.
9. **Prospective 3,804-patient evaluation grounding:** Incorporate the reported prospective evaluation in 3,804 oncology patients as support for oncology knowledge-graph grounding, neuro-symbolic multi-agent matching, criterion-level symbolic eligibility checks, explicit conflict resolution, audit trails, and clinician review loops; do not infer unstated benchmark or performance values.
10. **Prospective evaluation and operations:** Preserve knowledge-graph provenance, audit symbolic eligibility checks and criterion-level explanations, record multi-agent conflict resolutions and clinician review outcomes, and keep patient-level audit trails without inferring unreported performance benchmarks.
11. **Prospective validation pattern:** Use the reported 3,804-patient prospective evaluation as a validation design pattern: document criterion-level traceability, knowledge-graph provenance, agent-conflict resolution, clinician review, and prospective workflow evidence, while reporting only metrics supported by collected evidence.
12. **Primary implementation pattern:** Follow the prospective 3,804-patient evaluation pattern by grounding retrieval in an oncology-specific knowledge graph, assigning specialized agent roles, applying symbolic eligibility checks with explicit conflict resolution, attaching criterion-level evidence, validating prospectively, and producing clinician-auditable outputs.
13. **Reference architecture and validation:** Use the prospective 3,804-patient evaluation as a reference pattern for oncology knowledge-graph grounding, specialized agents, symbolic eligibility checks, criterion-level explanations, conflict resolution, and clinician review; evaluate local implementations only with locally collected evidence and without assuming unreported benchmark values.
14. **Prospective evidence scenario:** Treat the prospective 3,804-patient oncology trial-matching evaluation as an evidence scenario for neuro-symbolic eligibility reasoning, oncology knowledge-graph grounding, multi-agent role separation, conflict resolution, and clinician-auditable evidence outputs without inventing unreported performance values.
15. **Prospective 3,804-patient workflow:** Apply the reported workflow pattern by retrieving candidate trials through an oncology-specific knowledge graph, parsing eligibility with neuro-symbolic rules, resolving multi-agent conflicts explicitly, attaching criterion-level evidence, ranking patient-trial candidates for clinician review, and preserving clinician-auditable review loops.
16. **Validation anchor:** Use the published prospective 3,804-patient evaluation as a validation anchor for neuro-symbolic multi-agent oncology trial matching, documenting knowledge-graph-backed eligibility parsing, symbolic rule execution, conflict-resolution outcomes, patient-level audit trails, and prospective workflow reporting without adding unsupported benchmark claims.
17. **Prospective cohort reporting pattern:** In prospective matching runs, pair knowledge-graph-backed criterion parsing with separated specialist-agent roles, explicit conflict disposition, patient-level audit traces, and cohort-scale performance reporting limited to locally measured, source-supported evidence.
18. **Prospective performance reporting:** For local prospective deployments modeled on the 3,804-patient evaluation, report cohort size, data period, trial-source versions, knowledge-graph-backed criterion decisions, conflict-resolution dispositions, clinician-review outcomes, and measured performance only when actually collected; preserve criterion-level audit trails for every patient-trial decision.
19. **Prospective throughput and audit reporting:** In neuro-symbolic multi-agent matching deployments, ground retrieval in the oncology knowledge graph, preserve criterion-level reasoning and conflict-resolution dispositions, track cohort-scale throughput metrics across patient volumes such as the reported 3,804-patient evaluation, and provide clinician-auditable explanations for each patient-trial decision without inventing unreported performance values.
20. **Patient-scale clinician-auditable matching:** When adapting the reported 3,804-patient prospective evaluation pattern, combine oncology knowledge-graph eligibility reasoning, multi-agent conflict resolution, criterion-level evidence, and patient-scale performance reporting so clinicians can audit each output without relying on unsupported benchmark claims.
21. **Reference workflow for audit-ready matching:** Use the prospective 3,804-patient oncology trial-matching evaluation as a reference workflow that links oncology knowledge-graph retrieval, neuro-symbolic eligibility reasoning, multi-agent review, criterion-level evidence, explicit conflict resolution, and audit-ready clinician handoff without inferring unstated benchmark values.
22. **Prospective evidence integration:** Apply the 2026 prospective evidence from 3,804 oncology patients as support for neuro-symbolic multi-agent reasoning over an oncology-specific knowledge graph, requiring criterion-level evidence, explicit conflict-resolution records, scalability checks for thousands of patients, and clinician-auditable trial recommendations.
23. **Prospective evaluation practice package:** In prospective runs modeled on the 3,804-patient study, ground eligibility in an oncology-specific knowledge graph, retain criterion-level reasoning for each inclusion and exclusion decision, document multi-agent conflict resolution, preserve patient-trial audit trails, and report prospective performance metrics only when measured or source-supported.
24. **Benchmark for KG-backed eligibility reasoning:** Use the 2026 prospective 3,804-patient neuro-symbolic multi-agent oncology trial-matching evaluation as a benchmark pattern for knowledge-graph-backed eligibility reasoning, conflict resolution, criterion-level evidence, and clinician-auditable matching outputs without inventing unreported performance metrics.

## Inputs / Outputs

### Inputs

- Patient oncology summary or structured record with dates and source references.
- Current trial records containing identifiers, recruiting status, sites, interventions, and full eligibility criteria.
- Terminology mappings or an oncology knowledge graph for diagnoses, biomarkers, drugs, and procedures.
- Operational constraints such as travel radius, age range, trial phase, and preferred institutions.
- Optional clinician-provided interpretations or rules for ambiguous eligibility language.

### Outputs

- A ranked candidate-trial table with trial identifier, title, site, recruiting status, and match rationale.
- A criterion-level eligibility matrix labeled met, not met, unknown, or not applicable.
- Evidence links from each decision to patient facts and source trial text.
- A list of missing or stale data that could change eligibility.
- An agent-conflict log with the resolution, rationale, and items requiring clinician adjudication.
- A final status for each trial: likely eligible, potentially eligible pending data, likely ineligible, or insufficient evidence.
- An auditable run summary recording data versions, retrieval time, rules applied, and reviewer sign-off status.

Do not contact trial sites, alter patient care, or represent eligibility as confirmed. Verify current recruiting status and protocol details from authoritative trial sources, and require qualified clinician or trial-site confirmation before enrollment decisions.

## References

- Loaiza-Bonilla A, Yost C, Kurnaz S, Tuysuz E, Thaker NG. “Transforming oncology clinical trial matching through neuro-symbolic, multi-agent AI and an oncology-specific knowledge graph: a prospective evaluation in 3804 patients.” *ESMO Real World Data and Digital Oncology*. 2026 Jun. PubMed: https://pubmed.ncbi.nlm.nih.gov/42004487/
