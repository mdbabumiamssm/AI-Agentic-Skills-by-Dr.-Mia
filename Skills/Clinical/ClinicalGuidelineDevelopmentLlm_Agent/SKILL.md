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
7. **Live guideline-development evaluation design**: Retrieve evidence during drafting, preserve recommendation-level traceability, require expert adjudication, log disagreements and their disposition, measure evidence-retrieval, drafting, and review latency, and stop evaluation of any draft that is unsafe or unsupported.
8. **Worked prospective validation pattern**: Prospectively define and log each LLM-assisted task, its evidence inputs, outputs, elapsed time, and human actions; require expert adjudication against traceable evidence; grade each recommendation with the selected framework; calculate observed time savings against the established workflow without assuming a benchmark; analyze agreements, disagreements, and their dispositions; and require documented governance checkpoints before guideline adoption.
9. **Worked real-time governance pattern**: For guideline panels, capture each prompt, model/version, evidence input, draft output, recommendation grade, expert adjudication decision, and audit checkpoint so every recommendation can be traced from source evidence through LLM assistance to final panel disposition.
10. **Continuous recommendation lifecycle**: During active development, perform evidence surveillance; capture each source and retrieval time; draft or update recommendations only against captured evidence; log discrepancies between drafts, evidence, and current guidance with status and disposition; route unresolved or clinically material discrepancies to expert adjudication; preserve provenance for searches, evidence, LLM outputs, human edits, decisions, and approvals; and revise a recommendation when new evidence changes its certainty, benefit-harm balance, or applicability, when adjudication requires a change, when a logged discrepancy is upheld, or when its evidence or provenance becomes incomplete or outdated.
11. **Real-time guideline evaluation guardrails**: Use LLMs only for evidence organization, draft comparison, and consistency checks; log provenance for each recommendation; require expert adjudication for every graded recommendation; and screen failure modes including outdated evidence, omitted harms, and overconfident wording.
12. **Specialty-guideline concordance testing**: Compare LLM recommendations criterion by criterion against applicable professional-society recommendations, classifying discrepancies as omissions, contradictions, outdated advice, or potential harm.
13. **Conflict and dissent handling**: Track conflicts of interest, competing interpretations, panel disagreements, and the final governance disposition.
14. **Safety and scope control**: Flag unsupported claims, outdated evidence, population mismatches, off-label implications, and recommendations requiring human adjudication.
15. **Audit-ready output**: Produce a concise record of prompts, evidence inputs, LLM-generated drafts, human edits, decisions, and remaining uncertainties.
16. **Concurrent LLM evaluation checkpoints**: Evaluate each LLM-assisted guideline-development step in real time against traceable evidence and recommendation-grading support; require expert-panel checkpoint review, conflict disposition, and explicit documentation of where LLM output is accepted, revised, or rejected before it enters the guideline record.
17. **Real-time guideline-development governance pattern**: During clinical practice guideline drafting, require evidence traceability, expert adjudication, recommendation grading, bias checks, and explicit boundaries that limit LLM assistance to support tasks until accountable human reviewers approve inclusion in the guideline record.
18. **Real-time LLM guideline-development evaluation pattern**: During clinical practice guideline drafting, run live evidence checks for draft outputs, keep each recommendation traceable to source evidence and grading constraints, require expert adjudication before acceptance, review bias and hallucination risks, and preserve audit logs of prompts, evidence inputs, model outputs, reviewer decisions, and final disposition.
19. **Real-time safety adjudication pattern**: For each LLM-assisted guideline recommendation, keep evidence traceability, route the graded recommendation through expert adjudication, preserve audit logs of prompts, evidence inputs, model outputs, reviewer decisions, grades, and dispositions, and apply stopping rules when the recommendation is unsafe, unsupported, untraceable, or outside the approved evidence boundary.
20. **Real-time guideline text acceptance gate**: Before any LLM-generated guideline text is accepted, verify evidence traceability, confirm recommendation-grading support, capture expert adjudication, log conflicts or disagreements with their disposition, and complete the required governance checkpoint.

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
- Kabir R, Braud SC, Hinson CS, Nazerali RS. Are large language models consistent with the ASPS and AAPS guidelines? A comparison of AI chatbot recommendations and plastic surgery clinical guidance. *J Plast Reconstr Aesthet Surg*. 2026 May. https://pubmed.ncbi.nlm.nih.gov/41985209/
