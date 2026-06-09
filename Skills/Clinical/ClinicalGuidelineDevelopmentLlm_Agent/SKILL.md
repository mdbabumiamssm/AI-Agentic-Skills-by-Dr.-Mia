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
name: 'clinical-guideline-development-llm'
description: 'LLM-assisted clinical practice guideline development workflow grounded in real-time evaluation, evidence traceability, expert review, recommendation grading, and human governance.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Clinical Guideline Development LLM

## Overview

Use this skill to support clinical practice guideline development with an LLM while preserving evidence traceability, expert review, recommendation grading, conflict handling, and human governance. The workflow is grounded in real-time evaluation of LLM use for guideline development and is intended to assist, not replace, accountable clinical experts and guideline panels.

## When to Use This Skill

- Drafting or revising clinical practice guideline sections with explicit evidence links.
- Translating clinical questions into structured PICO-style evidence review tasks.
- Summarizing evidence tables, certainty judgments, and rationale for recommendations.
- Preparing recommendation language for expert panel review and grading.
- Reconciling reviewer comments, conflicts of interest, dissenting opinions, or unresolved evidence gaps.
- Auditing guideline text for traceability from recommendation to evidence source.
- Creating a governance-ready record of LLM assistance, human decisions, and unresolved risks.

## Core Capabilities

1. **Clinical question framing**: Convert a guideline topic into scoped questions, populations, interventions, comparators, outcomes, and decision context.
2. **Evidence traceability**: Maintain citations, study identifiers, review sources, and explicit links between evidence summaries and guideline statements.
3. **Recommendation drafting**: Generate guideline-ready recommendation text with rationale, caveats, implementation considerations, and review status.
4. **Recommendation grading support**: Organize certainty of evidence, balance of benefits and harms, patient values, feasibility, equity, and resource considerations for panel judgment.
5. **Real-time expert review workflow**: Present draft outputs in reviewable units, capture expert corrections, and revise only within approved evidence boundaries.
6. **Real-time evaluation protocol**: Require concurrent expert review of LLM outputs, trace each statement to its supporting evidence, classify detected errors, grade recommendations, adjudicate disagreements or uncertain content, and prevent any unreviewed LLM-generated text from entering the guideline.
7. **Conflict and dissent handling**: Track conflicts of interest, competing interpretations, panel disagreements, and the final governance disposition.
8. **Safety and scope control**: Flag unsupported claims, outdated evidence, population mismatches, off-label implications, and recommendations requiring human adjudication.
9. **Audit-ready output**: Produce a concise record of prompts, evidence inputs, LLM-generated drafts, human edits, decisions, and remaining uncertainties.

## Inputs / Outputs

**Inputs**

- Guideline topic, clinical scope, intended users, and target patient population.
- Clinical questions or draft PICO elements.
- Evidence sources, systematic reviews, study summaries, evidence tables, or PubMed/arXiv/GitHub links supplied by the user.
- Target grading framework, panel process, conflict-of-interest policy, and institutional governance requirements.
- Existing draft recommendations, reviewer comments, dissent statements, or implementation constraints.

**Outputs**

- Structured clinical question set and evidence map.
- Traceable evidence summary with source identifiers and unsupported-claim flags.
- Draft recommendations with rationale, grading placeholders, caveats, and review status.
- Reviewer-response table separating accepted edits, rejected edits, unresolved issues, and required expert decisions.
- Governance log documenting LLM assistance, human review, conflicts, dissent, and final disposition.

## References

- Erstad BL. Real-Time Evaluation of a Large Language Model for Clinical Practice Guideline Development. *Crit Care Explor*. 2026 May 1. https://pubmed.ncbi.nlm.nih.gov/42042855/
