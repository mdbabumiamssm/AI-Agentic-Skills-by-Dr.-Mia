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
name: 'ophthalmology-llm-safety-evaluation'
description: 'Evaluate ophthalmology LLM answers to CME or board-style questions for correctness, omissions, harm risk, guideline consistency, and clinician review.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Ophthalmology LLM Safety Evaluation

## Overview

Use this skill to evaluate large language model responses to ophthalmology continuing medical education, board-style, or clinical knowledge questions. The workflow emphasizes correctness, clinically important omissions, risk of harm, guideline consistency, model comparison, and transparent clinician-review requirements. It is intended for safety evaluation and research support, not autonomous diagnosis or treatment.

## When to Use This Skill

- Assess LLM answers to ophthalmology CME, board-review, resident teaching, or specialty clinical questions.
- Compare responses from multiple models on the same ophthalmology prompt set.
- Grade whether an answer is correct, partially correct, incomplete, misleading, or unsafe.
- Identify omitted findings, management steps, contraindications, referral thresholds, or follow-up intervals that matter in eye care.
- Evaluate whether a response aligns with current ophthalmology guidance, trial evidence, or specialty consensus.
- Prepare clinician-review forms, adjudication rubrics, or blinded evaluation summaries for ophthalmology LLM studies.
- Flag responses that could delay urgent care, misdirect treatment, or understate sight-threatening disease.

## Core Capabilities

1. Define the evaluation unit: Specify the prompt, expected answer source, model response, clinical topic, intended audience, and whether the task is CME, board-style, patient-facing, or clinician-facing.
2. Score correctness: Classify each response as correct, partially correct, incorrect, unsupported, or not assessable using the supplied answer key and cited clinical sources.
3. Detect content omission: Identify missing elements that materially affect diagnosis, management, prognosis, referral, counseling, or patient safety.
4. Rate risk of harm: Label harm risk as none, low, moderate, or high, with brief rationale tied to likely clinical consequences.
5. Check guideline consistency: Compare recommendations against current ophthalmology guidelines, drug labels, emergency referral standards, or authoritative review sources when available.
6. Capture specialty harm modes: Look for ophthalmology-specific safety failures, including missed acute angle closure, retinal detachment symptoms, endophthalmitis, giant cell arteritis, orbital cellulitis, chemical injury, optic neuritis, medication toxicity, amblyopia timing, and diabetic or glaucoma follow-up errors.
7. Support model comparison: Use the same rubric, source set, and adjudication rules across models, with anonymized response labels when possible.
8. Evaluate Gemini 3 Pro and GPT-5 family board-style benchmarks: For ophthalmology board-style question comparisons, record the exact model name/version and run date, stratify questions by specialty topic when available, score correctness, omissions, and harm with the same rubric across models, and limit benchmark use to education, evaluation, and clinician-supervised review rather than autonomous care.
9. Require clinician adjudication: Route uncertain, high-stakes, or discrepant assessments to an ophthalmologist or qualified clinician reviewer before treating them as final.
10. Produce audit-ready output: Return structured tables with prompt ID, topic, model label, correctness, omissions, harm rating, evidence notes, reviewer status, and adjudication comments.

## Inputs / Outputs

Inputs:

- Ophthalmology question, vignette, CME item, or board-style prompt.
- Model response or set of model responses to evaluate.
- Reference answer key, cited source, guideline, review article, textbook excerpt, or clinician-provided ground truth.
- Optional metadata: model name/version, temperature, prompt template, date run, audience, subspecialty, and reviewer identity or role.

Outputs:

- Structured evaluation table or JSON-like record for each prompt-response pair.
- Correctness label with concise rationale.
- Clinically important omission list.
- Risk-of-harm label with specialty-specific harm rationale.
- Guideline or evidence consistency notes with source links when available.
- Reviewer disposition: pending clinician review, reviewed, adjudicated, excluded, or unresolved.
- Summary of recurring error patterns across prompts or models, without inventing performance statistics.

## References

- Chen JL, Lu AJ, Verma R, Wang L, Koch DD. "Assessment of Correctness, Content Omission, and Risk of Harm in Large Language Model Responses to Ophthalmology Continuing Medical Education Questions." PubMed: https://pubmed.ncbi.nlm.nih.gov/41908501/
- Shean RS, Mallapu JK, Shah T, Rasheed HA, Younessi DN. "Comparative Performance of Gemini 3 Pro and GPT-5 Family Models on Ophthalmology Board-Style Questions." PubMed: https://pubmed.ncbi.nlm.nih.gov/41970036/
