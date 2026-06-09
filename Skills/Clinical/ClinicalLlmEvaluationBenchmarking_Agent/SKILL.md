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
   For domain-specific medical Q&A benchmarks, apply an item-level error taxonomy; check ambiguity, answerability, answer validity, clinical risk, and domain coverage; evaluate clinically relevant subgroup slices; verify answer alignment with supporting evidence; define and report reliability thresholds; and document dataset provenance, annotation decisions, omissions, hallucinations, leakage checks, prompt sensitivity, limitations, and trustworthy-use boundaries.

4. **Systematic clinical-note evaluation methods**  
   For AI-generated clinical note benchmarks, combine factual consistency checks, omission detection, risk-of-harm labels, note completeness scoring, and inter-rater reliability reporting; design benchmark datasets so automated evaluation methods can be experimentally compared with clinician scoring rather than assumed equivalent.

5. **Clinical-note quality measure taxonomy and validity**  
   Organize clinical-note quality measures into factuality, completeness, structure, readability, harmfulness, and clinician preference; report automated metrics separately from blinded human review, and assess benchmark validity by matching each measure to its intended quality dimension and experimentally checking whether automated results align with blinded clinical judgments.

6. **Metric selection**  
   Match the evaluation target to metric families such as factual consistency, clinical completeness, omission detection, guideline concordance, readability, usefulness, calibration, and safety risk. Treat lexical similarity and embedding similarity as supporting signals rather than standalone clinical quality measures.

7. **Clinician adjudication design**  
   Define reviewer qualifications, rating rubrics, blinded review procedures, disagreement resolution, inter-rater agreement reporting, escalation criteria, and examples for each label.

8. **Omission and harm scoring**  
   Score missing critical information separately from incorrect additions, and classify likely harm severity using task-specific clinical risk categories.

9. **Safety and bias labeling**  
   Label hallucinations, unsupported claims, contradictions, unsafe triage, medication errors, delayed-care risk, inappropriate certainty, protected-class bias, and patient comprehension risks.

10. **Experimental benchmark execution**  
   Compare candidate systems under controlled prompts, fixed input sets, reproducible model settings, documented retrieval sources, and consistent post-processing.

11. **Regression gates and release criteria**  
   Convert evaluation results into pass/fail thresholds for deployment, including minimum clinician approval, maximum critical-error rate, no unresolved severe harm findings, and no regression on sentinel cases.

12. **Result interpretation**  
   Summarize quantitative scores with confidence intervals where appropriate, describe representative failures, separate clinical risk from style preference, and document limits of generalization.

13. **Evaluation artifact packaging**  
   Produce a reusable benchmark plan, annotation guide, data card, metric table, adjudication summary, risk register, regression checklist, and deployment recommendation.

14. **On-premises diagnosis deployment evaluation**  
   For open-source distilled reasoning LLMs proposed for on-premises clinical diagnosis, design controlled evaluations of diagnostic accuracy, calibration, latency, privacy, failure modes, quantization effects, and local hardware constraints against appropriate baselines; define task-specific deployment go/no-go thresholds and require safety review, governance approval, and regression gates before clinical use.

15. **ICU nursing QA evaluation profile**  
   Construct specialty-specific intensive care nursing question sets and run blinded comparisons across proprietary and open reasoning models; evaluate factual accuracy, clinical relevance, and clinically important omissions separately, and define escalation requirements for uncertain, unsafe, or high-risk answers.

16. **Clinical response correctness and harm evaluation**  
   Evaluate correctness, clinically important omission, and risk of harm as separate dimensions using clinician adjudication, severity-weighted scoring, and confidence intervals; report explicit examples of unsafe answers rather than relying on aggregate accuracy alone.

17. **Fine-grained domain Q&A adjudication**  
   Evaluate domain-specific clinical Q&A at the item level using clinically meaningful error taxonomies, separate omission and harm labels, clinically relevant subgroup analysis, evidence-grounding checks, and calibration assessment; adjudicate individual answers and failure modes rather than relying on aggregate accuracy alone.

18. **Trustworthy domain-specific Q&A dataset auditing**  
   Audit each benchmark item for ambiguity, evidence sufficiency, omissions, harmfulness, subgroup difficulty, and potential contamination; report these findings transparently alongside benchmark limitations and trustworthy-use boundaries.

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
- PubMed: Gülhan Güner S, Tan Z, Gülpınar S. "Comparative performance of artificial intelligence models in intensive care nursing questions: an evaluation of ChatGPT, DeepSeek, and Google Gemini." *BMC Nursing*. 2026 May 2. https://pubmed.ncbi.nlm.nih.gov/42069581/
- PubMed: Chen JL, Lu AJ, Verma R, Wang L, Koch DD. "Assessment of Correctness, Content Omission, and Risk of Harm in Large Language Model Responses to Ophthalmology Continuing Medical Education Questions." *Ophthalmol Sci*. 2026 May. https://pubmed.ncbi.nlm.nih.gov/41908501/
