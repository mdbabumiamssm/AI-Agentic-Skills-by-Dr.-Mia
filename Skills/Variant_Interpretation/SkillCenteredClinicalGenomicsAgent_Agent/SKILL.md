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
name: 'skill-centered-clinical-genomics-agent'
description: 'Orchestrate a research workflow for genomic interpretation using task-specific skills, traceable guideline and database evidence, computational analysis, and specialist review, informed by ClinGeX.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Skill-Centered Clinical Genomics Agent

## Overview

This skill organizes research genomic interpretation into selected, dependency-aware operations with traceable evidence and explicit expert review. It is informed by ClinGeX, which combines 105 curated executable skills through task-specific retrieval and orchestration to address fragmented genomic evidence. This document defines an adaptation for research evaluation; it neither installs ClinGeX nor establishes validated clinical performance.

## When to Use This Skill

- Assemble evidence for criterion-level variant interpretation or rare-disease genomic research.
- Investigate gene–disease validity or molecular oncology evidence in a defined clinical context.
- Coordinate existing guideline, retrieval, annotation, and prediction tools with explicit prerequisites.
- Evaluate concordance with frozen expert reference classifications and inspect evidence provenance.
- Prepare an auditable evidence packet for a qualified clinical genomics specialist.

## Core Capabilities

1. **Extract the source architecture before implementation.** Consult the versioned paper, especially Methods 4.1 and Supplementary A.3–A.4, and record the sections used. The reported library contains 27 guide, 40 search, and 38 analysis skills. Its selector decomposes each message into one to three needs, ranks skills by maximum embedding similarity, selects 30 seeds, adds up to 15 one-hop related skills with relationship score at least 8, includes three fixed literature-search skills, and closes hard prerequisites. Instructions are progressively disclosed; execution records preserve actions and outputs. Treat these as paper-specific design choices, not settings validated for a different library. Do not claim reproduction unless the corresponding artifacts and configuration are available.

2. **Separate reported evidence from local success criteria.** The supplied finding reports 71.7% macro-averaged task-level concordance across 223 temporally held-out cases, versus 45.7% for the strongest comparator. Its counterbalanced collaboration study involved 45 professionals and 135 cases; reported ClinGeX concordance was 51.1%, 57.8%, and 64.4% across the three clinical scenarios. The full text identifies the abstract's corresponding 15.6%, 31.1%, and 20.0% comparator values as DeepSeek-assisted results. The exploratory month involved 11 professionals and 96 self-selected episodes. These are study results, not expected accuracy for this skill, diagnostic accuracy in an unselected population, or evidence of improved patient outcomes. Preserve distinctions between automated task-level metrics, human case-level concordance, and self-reported utility. See the source finding and full text in References.

3. **Define a bounded research question.** The following steps are proposed operating instructions for this adaptation. Select one task and document its scope, evidence cutoff, and intended output. Confirm genome assembly, chromosome and alleles, transcript accession/version, variant representation, and disease context before interpreting a variant. Record phenotype, inheritance, segregation, specimen type, and tumor context when relevant. Request missing indispensable identifiers or return an explicitly incomplete result; never silently choose an assembly or transcript. Use de-identified inputs and avoid transmitting identifiable clinical records to public retrieval services.

4. **Select available skills and resolve dependencies.** Inspect the actual local inventory; empty inventories do not imply that the paper's skills are installed. For each selected operation, record its identifier/version, required inputs, prerequisites, executable interface, and expected output. Use semantic retrieval when an index exists; otherwise document manual task-to-capability selection. Build a directed prerequisite plan, detect cycles, and execute prerequisites before dependents. Report unavailable skills as gaps and do not fabricate tool results or substitute an undocumented algorithm.

5. **Retrieve guidelines and database evidence with provenance.** Identify the applicable current guideline and any gene- or disease-specific specification, recording version and retrieval date. Retrieve relevant primary literature and database records, such as ClinVar, ClinGen, population frequency resources, or oncology knowledge bases, according to the task. Record accession, URL, release/version, relevant passage or data field, applicability, review status, and disagreements. Distinguish a database assertion from the underlying evidence; multiple records citing one study are not independent support. Unreachable sources remain unavailable, rather than negative evidence.

6. **Integrate computational evidence without overcounting.** Run only available, inspected operations on normalized inputs. Preserve software/model version, parameters, reference assembly, input identifiers, output, warnings, and execution status. Assess applicability of population, consequence, splice, functional, phenotype, and segregation evidence separately. Apply only thresholds supported by the selected guideline or validated model specification. Do not count correlated predictions as independent evidence, infer pathogenicity from a score alone, or equate a predicted molecular effect with clinical actionability.

7. **Synthesize an auditable draft.** Map each candidate criterion to evidence identifiers, the applicable rule, evidence strength, rationale, and a status such as supported, not met, conflicting, or not assessable. Combine criteria only under the selected framework. Keep germline pathogenicity, gene–disease validity, and tumor-context actionability as distinct outputs; do not translate their labels interchangeably. For oncology, retain tumor type, indication, treatment setting, evidence level, and evidence date. If essential dependencies fail, limit conclusions to supported observations and list the unresolved questions.

8. **Review and evaluate reproducibly.** Route the draft, conflicting evidence, and unsupported assumptions to a qualified specialist, retaining pending review status until review occurs. Do not issue an autonomous diagnosis, treatment instruction, or clinical sign-out. For benchmark use, freeze case inputs, source versions, reference labels, and scoring rules before execution; prevent direct or indirect access to held-out answer records while permitting the prespecified underlying evidence. Log retrievals to audit leakage. Report per-task denominators, concordance, abstentions, failed runs, and elapsed time; state macro versus pooled aggregation explicitly. Keep benchmark retrieval restrictions distinct from ordinary evidence synthesis, where curated assertions may be legitimate inputs.

## Inputs / Outputs

**Inputs**

- A de-identified question and task type: variant criteria, rare disease, gene–disease validity, or oncology evidence.
- Normalized variant identifiers and assembly/transcript versions, plus relevant phenotype, inheritance, family, specimen, or tumor context.
- An available skill/tool inventory, accessible evidence sources, and guideline versions or instructions to retrieve them.
- For evaluation: a frozen case set, evidence cutoff, independently held reference answers, and prespecified scoring rules.

**Outputs**

Return a structured research packet using these fields:

| Field | Required contents |
| --- | --- |
| `run` | Task, de-identified case ID, scope, start/end times, elapsed time, evidence cutoff, completion status |
| `source_basis` | Paper version/sections, extracted architecture, reported evaluation, departures from the source design |
| `normalized_inputs` | Assembly, transcript, variant representation, clinical context, missing inputs |
| `execution_plan` | Selected skill versions, dependencies, execution order, status, failed or unavailable operations |
| `evidence_ledger` | Evidence ID, source/accession, URL, release/retrieval date, finding, applicability, dependencies and conflicts |
| `computational_results` | Tool/model versions, parameters, inputs, outputs, warnings and limitations |
| `criterion_matrix` | Framework/version, criterion, strength, evidence IDs, rationale and assessment status |
| `draft_interpretation` | Supported research conclusion, alternatives, uncertainty, unresolved questions and scope limits |
| `expert_review` | Pending/reviewed status, specialist decisions and rationale for changes |
| `evaluation` | If applicable: scoring rule, task denominators, concordance, failures, abstentions and leakage audit |

**Completion check and time budget**

Aim to produce the packet within 15 minutes for a bounded case with available dependencies. Check that every conclusion resolves to evidence and an applicable rule, all essential identifiers are explicit, computational claims have actual execution records, and review status is accurate. The time target is an operational objective, not a result established by the paper or a deadline for clinical sign-out. If retrieval or execution prevents completion, return a partial or blocked packet within the budget, specify the next required action, and do not count it as successful interpretation.

## References

- [Advancing Clinical Genomic Analysis with a Skill-Centred Agent — source finding, arXiv:2610.40556v1](http://arxiv.org/abs/2610.40556v1). Research preprint; source of the reported study outcomes.
- [Versioned full text and supplementary methods](https://arxiv.org/html/2610.40556v1). Sections 2.1–2.4, 4.1, A.3, and A.4 provide architecture, retrieval, evaluation, and scoring context. Clinical evidence and guideline versions must be retrieved separately for each research task.
