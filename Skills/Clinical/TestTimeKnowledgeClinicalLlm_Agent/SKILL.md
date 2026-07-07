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
name: 'test-time-knowledge-clinical-llm'
description: 'Use test-time knowledge acquisition to retrieve, vet, inject, and cite current clinical evidence for LLM-assisted medical decision support without fine-tuning.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Test-Time Knowledge Clinical LLM

## Overview

This skill guides clinical LLM workflows that acquire and integrate medical knowledge during inference rather than relying only on static training data or fine-tuning. It emphasizes trustworthy retrieval, source ranking, citation-grounded reasoning, conflict handling, safety checks, and explicit clinician oversight for medical decision support.

Use this skill to make a clinical answer traceable to current evidence while preserving uncertainty, scope limits, and human accountability.

## When to Use This Skill

- A clinical LLM response needs current medical evidence retrieved at test time.
- The task involves medical decision support, differential diagnosis, triage reasoning, therapy comparison, guideline lookup, or patient-specific evidence synthesis.
- The model must cite sources and distinguish retrieved evidence from internal model knowledge.
- Evidence may conflict across guidelines, studies, drug labels, clinical reviews, or institutional protocols.
- The answer could affect diagnosis, treatment, medication dosing, contraindication screening, escalation, or patient safety.
- The user asks for retrieval-augmented, evidence-grounded, test-time adaptation, or knowledge-base-assisted clinical reasoning.

## Core Capabilities

1. **Clinical query framing**: Convert the user request into a focused clinical question, including population, condition, intervention, comparator, outcome, urgency, and setting when available.

2. **Retrieval scope control**: Select evidence sources appropriate to the decision, prioritizing current guidelines, regulatory labels, systematic reviews, clinical trials, high-quality reviews, and institution-approved material.

3. **Source trust ranking**: Rank sources by clinical authority, recency, applicability, evidence type, transparency, and patient-context fit before using them in reasoning.

4. **Test-time knowledge injection**: Insert concise retrieved facts, guideline statements, contraindications, uncertainty notes, and source citations into the model context before generating a clinical response.

5. **Conflict handling**: Identify disagreements across sources, explain likely reasons such as population differences or publication date, and avoid forcing consensus when evidence remains uncertain.

6. **Citation grounding**: Tie clinically important claims to citations, especially recommendations, contraindications, diagnostic thresholds, drug facts, and risk statements.

7. **Safety and escalation checks**: Screen for emergency symptoms, high-risk medication issues, pregnancy, pediatrics, renal or hepatic impairment, immunocompromise, allergy, and other factors requiring clinician review.

8. **Uncertainty management**: Separate known evidence, inferred reasoning, missing information, and questions that require additional clinical data or direct evaluation.

9. **Clinician oversight workflow**: Present outputs as decision support, not autonomous care, and require qualified clinical review before diagnosis, treatment, or workflow action.

10. **Knowledge-base hygiene**: When maintaining a reusable test-time knowledge base, add, update, merge, or retire entries with provenance, date checked, clinical scope, and applicability limits.

11. **Test-time acquisition procedure and evaluation**: For each clinical question, generate focused retrieval subqueries and multiple reasoning trajectories; select only decision-relevant evidence after checking source authority, recency, applicability, and provenance; synthesize a pseudo-reference answer and use self-consistency to separate confident cases for reward-guided heuristic extraction from unconfident cases for reflection on reasoning gaps. Reconcile contradictions explicitly, bind every material claim and acquired knowledge entry to its supporting citation, and abstain or escalate when evidence is insufficient, irreconcilable, or unsafe. Add, modify, or merge validated knowledge in a capacity-controlled knowledge base that guides later inference without parameter updates, and evaluate the workflow against both zero-shot and static-RAG baselines without asserting improvements unless measured.

12. **Acquisition controls and ablation**: Control test-time query generation for clinical specificity and coverage; rank source authority before inclusion; enforce evidence-freshness checks against the decision date; preserve, reconcile, or escalate material contradictions; allocate the context budget by decision relevance, authority, and safety impact; verify that each citation entails the claim it supports; and ablate retrieved knowledge to identify conclusions that depend on unsupported parametric reasoning rather than acquired evidence.

13. **Test-time knowledge-acquisition runbook**: Generate targeted queries from the clinical question and unresolved reasoning gaps; check publication, revision, and access dates for freshness; vet source authority, provenance, applicability, and evidentiary role before inclusion; surface conflicts and resolve them only when source differences justify doing so; stop retrieval when decision-critical claims are adequately supported, new searches add no material evidence, or time, context, safety, or source limits require abstention or escalation; bind each material claim to the evidence that supports it; and ablate the workflow against a no-retrieval baseline to identify which outputs depend on acquired knowledge without claiming improvement unless measured.

14. **Evidence-controlled acquisition pipeline**: Generate focused queries for the clinical decision and identified knowledge gaps; vet retrieved sources for authority, provenance, recency, and patient-context applicability; deduplicate overlapping evidence by claim and source lineage; construct a bounded context that preserves decision-critical evidence, limitations, and citations; handle conflicts by reporting unresolved disagreement or preferring a source only when its authority, currency, or applicability justifies that choice; verify that citations resolve to the intended source and support the associated claim; run acquisition ablations such as no-retrieval, source-removal, and evidence-removal checks without asserting performance gains unless measured; enforce search-count, time, and context-size limits to control latency; and fail closed by withholding the affected recommendation and escalating for clinician review when trustworthy evidence is unavailable, unverifiable, materially conflicting, or insufficient.

15. **Inference-time evidence acquisition guardrails**: Retrieve current evidence at inference time for each medical decision; grade and filter candidate sources before injecting them into context; cite the specific evidence spans supporting material claims; label static model knowledge separately from retrieved knowledge; and evaluate whether retrieval improves correctness without adding unsafe authority bias or overconfident deference to retrieved sources.

16. **Retrieval-and-verification deployment pattern**: For medical decision-making tasks, generate focused test-time queries, vet candidate sources before use, inject only verified evidence into the response context, capture citations for every material claim, check evidence freshness against the clinical decision date, and keep the workflow within no-fine-tune deployment boundaries so acquired knowledge informs inference without changing model parameters.

17. **Decision-triggered knowledge acquisition**: Trigger retrieval when a medical decision support answer depends on current guidelines, drug labels, diagnostic thresholds, safety warnings, patient-specific applicability, or uncertain model knowledge; vet sources for authority, provenance, publication or revision date, applicability, and evidentiary role before injection; format injected evidence as claim, source, date checked, applicability, limitations or contradictions, and citation; cite every material recommendation, contraindication, threshold, drug fact, and risk statement; compare source dates against the clinical decision date; preserve unresolved contradictions instead of forcing consensus; and abstain or escalate when retrieved evidence is weak, unverifiable, outdated, materially conflicting, or insufficient for the requested decision.

18. **Unsupported-reasoning safeguards**: For clinical decision support, start test-time acquisition whenever recommendations, contraindications, diagnostic thresholds, drug facts, patient-specific applicability, or model uncertainty require current evidence; vet each retrieved source for authority, provenance, freshness, applicability, and direct support before citation injection; include dates checked and applicability limits with injected evidence; and block, qualify, or escalate any conclusion that remains uncited, stale, unverifiable, or only supported by internal model reasoning.

19. **Test-time retrieval-and-verification workflow**: Plan retrieval queries from the clinical decision and unresolved knowledge gaps; vet candidate sources for authority, provenance, recency, applicability, and evidentiary role; inject only verified evidence with captured citations and freshness dates into the response context; and guard against ungrounded reasoning by blocking, qualifying, or escalating unsupported conclusions without fine-tuning or parameter updates.

## Inputs / Outputs

**Inputs**

- Clinical question, task, case vignette, or patient-specific context.
- Available patient factors such as age, sex, pregnancy status, comorbidities, medications, allergies, labs, imaging, vitals, and care setting.
- Retrieval constraints such as approved sources, jurisdiction, specialty, guideline body, publication date range, or institutional policy.
- Output constraints such as required format, citation style, risk tolerance, urgency, or clinician-facing versus patient-facing language.

**Outputs**

- A concise clinical answer grounded in retrieved evidence.
- Source-ranked evidence summary with citation links and dates when available.
- Clear distinction between retrieved facts, model reasoning, assumptions, and unresolved uncertainty.
- Conflict notes when credible sources disagree.
- Safety checks, red flags, contraindications, and escalation recommendations.
- Follow-up questions needed for safe interpretation.
- Clinician oversight statement for any diagnosis, treatment, medication, or operational recommendation.

## References

- Li S, Bao L, Li S, Wan B. Enhancing LLM-based medical decision-making by test-time knowledge acquisition. PubMed PMID: 41953846. https://pubmed.ncbi.nlm.nih.gov/41953846/
- Lewis P, Perez E, Piktus A, et al. Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks. arXiv:2005.11401. https://arxiv.org/abs/2005.11401
- SurgeryLLM: a retrieval-augmented generation large language model framework for surgical decision support and workflow enhancement. PubMed PMID: 39695316. https://pubmed.ncbi.nlm.nih.gov/39695316/
