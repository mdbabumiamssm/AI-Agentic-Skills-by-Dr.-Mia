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
