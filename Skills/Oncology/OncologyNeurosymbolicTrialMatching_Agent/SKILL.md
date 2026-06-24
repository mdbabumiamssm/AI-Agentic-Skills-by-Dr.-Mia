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
9. **Prospective evaluation and operations:** Preserve knowledge-graph provenance, audit symbolic eligibility checks and criterion-level explanations, record multi-agent conflict resolutions and clinician overrides, audit recall against reviewed relevant trials, and track operational throughput metrics without inferring unreported performance benchmarks.
10. **Prospective validation pattern:** Use the reported 3,804-patient prospective evaluation as a validation design pattern: measure criterion-level traceability, knowledge-graph provenance coverage, agent-conflict resolution outcomes, turnaround time, clinician overrides, and prospective matching metrics, while reporting only metrics supported by collected evidence.
11. **Primary implementation pattern:** Follow the prospective 3,804-patient evaluation pattern by grounding retrieval in an oncology-specific knowledge graph, assigning specialized agent roles, applying symbolic eligibility checks with explicit conflict resolution, attaching criterion-level evidence, tracking coverage metrics, validating prospectively, and producing clinician-auditable outputs.
12. **Reference architecture and validation:** Use the prospective 3,804-patient evaluation as a reference pattern for oncology knowledge-graph grounding, specialized agents, symbolic eligibility checks, criterion-level explanations, conflict resolution, and clinician review; evaluate recall and reviewer workload using locally collected evidence without assuming unreported benchmark values.
13. **Benchmark scenario:** Treat the prospective 3,804-patient oncology trial-matching evaluation as a benchmark scenario for neuro-symbolic eligibility reasoning, oncology knowledge-graph grounding, multi-agent role separation, conflict resolution, and clinician-auditable evidence outputs without inventing unreported performance values.

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
