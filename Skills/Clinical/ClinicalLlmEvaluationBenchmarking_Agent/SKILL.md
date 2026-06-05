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
name: 'clinical-llm-evaluation-benchmarking'
description: 'Design and run rigorous clinical LLM evaluation workflows grounded in systematic review evidence on AI-generated clinical note quality methods.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Clinical LLM Evaluation Benchmarking

## Overview

This skill guides the design, execution, and interpretation of evaluations for clinical LLM outputs, with emphasis on AI-generated clinical notes and adjacent clinical text workflows. It helps select appropriate metrics, clinician review procedures, safety labels, omission and harm scoring, dataset provenance checks, and regression gates before any clinical deployment or workflow expansion.

Use this skill to avoid relying on a single automatic score when clinical usefulness, factual correctness, patient safety, and documentation risk all need to be evaluated together.

## When to Use This Skill

- Evaluating AI-generated clinical notes, summaries, discharge instructions, consult drafts, radiology-style text, or handoff documentation.
- Comparing multiple LLMs, prompts, retrieval pipelines, note templates, or deployment configurations for clinical use.
- Designing a benchmark dataset with known provenance, consent constraints, de-identification status, and clinical specialty coverage.
- Creating safety labels for hallucination, omission, contradiction, bias, inappropriate advice, unsafe escalation, or guideline discordance.
- Setting clinician adjudication workflows for factuality, completeness, usefulness, readability, and potential patient harm.
- Building regression gates before model, prompt, retrieval corpus, or EHR integration changes.
- Auditing whether automatic metrics align with expert clinical judgment for a specific clinical task.

## Core Capabilities

1. **Evaluation scope definition**  
   Identify the clinical task, intended users, deployment context, acceptable use boundaries, and failure modes before choosing metrics.

2. **Dataset provenance and governance review**  
   Record source systems, cohort criteria, date ranges, de-identification method, annotation status, specialty mix, inclusion and exclusion criteria, and patient privacy constraints.

3. **Fine-grained medical Q&A dataset evaluation**  
   For domain-specific medical Q&A benchmarks, document question quality, answer validity, clinical risk, domain coverage, item-level error taxonomies, omission and hallucination labels, leakage checks, prompt sensitivity checks, and trustworthy reporting templates needed to support model comparisons.

4. **Systematic clinical-note evaluation methods**  
   For AI-generated clinical note benchmarks, combine factual consistency checks, omission detection, risk-of-harm labels, note completeness scoring, and inter-rater reliability reporting; design benchmark datasets so automated evaluation methods can be experimentally compared with clinician scoring rather than assumed equivalent.

5. **Metric selection**  
   Match the evaluation target to metric families such as factual consistency, clinical completeness, omission detection, guideline concordance, readability, usefulness, calibration, and safety risk. Treat lexical similarity and embedding similarity as supporting signals rather than standalone clinical quality measures.

6. **Clinician adjudication design**  
   Define reviewer qualifications, rating rubrics, blinded review procedures, disagreement resolution, inter-rater agreement reporting, escalation criteria, and examples for each label.

7. **Omission and harm scoring**  
   Score missing critical information separately from incorrect additions, and classify likely harm severity using task-specific clinical risk categories.

8. **Safety and bias labeling**  
   Label hallucinations, unsupported claims, contradictions, unsafe triage, medication errors, delayed-care risk, inappropriate certainty, protected-class bias, and patient comprehension risks.

9. **Experimental benchmark execution**  
   Compare candidate systems under controlled prompts, fixed input sets, reproducible model settings, documented retrieval sources, and consistent post-processing.

10. **Regression gates and release criteria**  
   Convert evaluation results into pass/fail thresholds for deployment, including minimum clinician approval, maximum critical-error rate, no unresolved severe harm findings, and no regression on sentinel cases.

11. **Result interpretation**  
   Summarize quantitative scores with confidence intervals where appropriate, describe representative failures, separate clinical risk from style preference, and document limits of generalization.

12. **Evaluation artifact packaging**  
   Produce a reusable benchmark plan, annotation guide, data card, metric table, adjudication summary, risk register, regression checklist, and deployment recommendation.

13. **On-premises diagnosis deployment evaluation**  
   For open-source and distilled reasoning LLMs proposed for clinical diagnosis, design controlled benchmarks that compare diagnostic performance, latency, privacy tradeoffs, local hardware and deployment constraints, calibration, and failure modes against appropriate baselines; require safety review, governance approval, and performance/regression gates before clinical use.

## Inputs / Outputs

**Inputs**

- Clinical LLM use case, target users, intended care setting, and non-use boundaries.
- Candidate model outputs, prompts, retrieval settings, system configuration, and comparison baselines.
- Benchmark dataset or dataset specification, including provenance and governance constraints.
- Reference notes, source records, guidelines, gold annotations, or clinician-authored comparison material when available.
- Safety taxonomy, rubric requirements, reviewer pool, and acceptable clinical risk thresholds.

**Outputs**

- Clinical LLM evaluation plan with task scope, hypotheses, metrics, and reviewer workflow.
- Dataset provenance and readiness checklist.
- Clinician adjudication rubric with label definitions and severity levels.
- Automatic and human evaluation matrix matched to the clinical task.
- Omission, hallucination, contradiction, and harm scoring summary.
- Regression gate checklist for model, prompt, retrieval, or integration changes.
- Final benchmark report with deployment recommendation, residual risks, and required follow-up validation.

## References

- PubMed: Dahlberg A, Käenniemi T, Winther-Jensen T, Tapiola O, Luisto R. "Measuring the quality of AI-generated clinical notes: A systematic review and experimental benchmark of evaluation methods." *Artificial Intelligence in Medicine*. 2026 Jul. https://pubmed.ncbi.nlm.nih.gov/41955894/
- PubMed: Zhong W, Fu Y, Peng D, Liu Y, Liu Y. "Open-Source Large Language Models Distilled DeepSeek-R1 Pose Challenges for On-Premises Clinical Deployment in Medical Diagnosis: A Comparative Study of Performance." *Journal of Medical Systems*. 2026 May 1. https://pubmed.ncbi.nlm.nih.gov/42062641/
- PubMed: Fonseca RDC, Rios RA, Castaldoni R, Carvalho AA, Lopes TJS. "Fine-grained evaluation of a domain-specific Q&A dataset to support trustworthy medical language models." *Health Inf Sci Syst*. 2026 Dec. https://pubmed.ncbi.nlm.nih.gov/42039929/
