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
   For open-source distilled reasoning LLMs proposed for on-premises clinical diagnosis, design controlled evaluations of diagnostic accuracy, calibration, uncertainty reporting, quantization and distillation degradation, local hardware constraints, privacy-constrained deployment checks, latency, and failure modes; select clinically relevant baselines, including proprietary-model comparators when appropriate; define deployment-specific safety gates and task-specific go/no-go thresholds versus those baselines; and require safety review, governance approval, and regression testing before clinical use.

15. **Specialty-stratified cross-provider ICU nursing QA benchmarking**  
   Construct specialty-specific intensive care nursing question sets and run blinded comparisons across providers such as ChatGPT, DeepSeek, and Google Gemini using expert-reference scoring; apply appropriate statistical model comparisons, report uncertainty, evaluate factual accuracy, clinical relevance, and clinically important omissions separately, and define escalation requirements for uncertain, unsafe, or high-risk answers. Treat exam-question performance as task-specific evidence and do not generalize it to bedside safety without prospective clinical workflow validation.

16. **Clinical response correctness and harm evaluation**  
   Evaluate correctness, clinically important omission, and risk of harm as separate dimensions using clinician adjudication, severity-weighted scoring, and confidence intervals; report explicit examples of unsafe answers rather than relying on aggregate accuracy alone.

17. **Fine-grained domain Q&A adjudication**  
   Evaluate domain-specific clinical Q&A at the item level using clinically meaningful error taxonomies, separate omission and harm labels, clinically relevant subgroup analysis, evidence-grounding checks, and calibration assessment; adjudicate individual answers and failure modes rather than relying on aggregate accuracy alone.

18. **Trustworthy domain-specific Q&A dataset auditing**  
   Audit each benchmark item for ambiguity, evidence sufficiency, omissions, harmfulness, subgroup difficulty, and potential contamination; report these findings transparently alongside benchmark limitations and trustworthy-use boundaries.

19. **Fine-grained domain Q&A audit controls**  
   Audit domain-specific Q&A datasets with an item taxonomy, ambiguity and answerability checks, clinician adjudication, separate omission and harmfulness labels, subgroup analysis, and controls against benchmark contamination.

20. **Blinded specialist adjudication of correctness, omission, and harm**  
   Evaluate correctness, clinically important omission, and risk of harm as separate endpoints; require blinded specialist adjudication with a prespecified process for resolving disagreements; weight errors by clinical severity; and separately analyze unsafe answers stated with confidence.

21. **Atomic scoring for fine-grained domain Q&A evaluation**
   Score domain-specific Q&A outputs at the atomic item or answer-component level for correctness, completeness, evidence support, uncertainty, harmful omissions, and subgroup performance; audit dataset quality alongside model performance, and use aggregate accuracy only as a secondary summary rather than the sole evaluation criterion.

22. **Systematic clinical-note quality measurement and reporting**
   Separate clinical-note quality measures into factuality, omission, structure, usability, safety, and downstream task utility; report evaluator reliability, assess agreement between human judgments and automated metrics, and use benchmark reporting templates that document each measure, evaluator, comparison method, and result.

23. **Fine-grained medical Q&A evaluation pattern**
   Build domain-specific Q&A datasets with documented construction criteria; label omissions and hallucinations separately; apply an item-level clinical error taxonomy; report trustworthiness metrics for answer validity, evidence support, clinical risk, and reliability; and use benchmark reporting templates that disclose dataset provenance, annotation decisions, limitations, and trustworthy-use boundaries.

24. **Clinical education benchmark module**
   For board-style or CME question benchmarks, score correctness, content omission, and risk of harm as independent endpoints; prespecify clinician adjudication rules for reviewer qualifications, disagreement resolution, and severity handling; and report results by specialty rather than relying only on pooled aggregate performance.

25. **Fine-grained domain-specific medical Q&A dataset controls**
   For trustworthy medical LLM benchmarking, evaluate domain-specific Q&A datasets with item-level error taxonomies, separate omission and harm labels, domain coverage checks, calibration by question type, and dataset provenance controls.

26. **Fine-grained domain Q&A evaluation workflow**
   Construct domain-specific Q&A datasets with documented source, inclusion, and evidence criteria; use clinician adjudication to label ambiguity, omissions, evidence grounding, and item-level error categories; assess calibration and clinically relevant subgroups; and report item-level findings alongside aggregate accuracy.

27. **Intensive-care nursing comparative evaluation**
   For intensive-care nursing questions, use expert reference answers, direct model-to-model comparison, clinically weighted error categories, and dedicated analysis of unsafe omissions; report these dimensions separately rather than relying on aggregate accuracy alone.

28. **Fine-grained medical Q&A dataset quality controls**
   Label item-level difficulty and ambiguity, assess answer-evidence alignment, apply separate omission and contradiction taxonomies, report clinically relevant subgroup results and annotator agreement, and check for dataset contamination.

29. **Management-critical omission and harm adjudication**
   Score correctness, content omission, and risk of harm as separate endpoints; require severity-weighted clinician adjudication; and explicitly report cases where fluent answers omit management-critical details.

30. **Fine-grained domain Q&A evidence and reliability audit**
   Evaluate each item with a prespecified error taxonomy; score omissions separately; verify supporting evidence; flag ambiguity and define handling rules; report clinically relevant subgroup results and inter-rater agreement; and check benchmark items and source data for dataset leakage.

31. **Specialty-aware intensive-care nursing benchmark interpretation**
   For intensive-care nursing questions, use clinician-authored scoring rubrics and reference answers, compare ChatGPT, DeepSeek, and Google Gemini under the same question set and conditions, and analyze safety-critical errors and omissions separately from overall accuracy. Treat board-style question accuracy as task-specific evidence that cannot establish bedside performance or safety without validation in clinical workflows.

32. **Fine-grained medical Q&A validation and failure analysis**
   Validate domain-specific medical Q&A datasets with an item-level taxonomy; score omissions and uncertainty separately; analyze clinically relevant subgroups; and explicitly check whether aggregate accuracy conceals clinically important failures.

33. **Ophthalmology CME answer safety adjudication**
   Score correctness, clinically material omission, and risk of harm separately; require blinded expert adjudication; and report harmful-error rates independently from overall answer accuracy.

34. **Versioned fine-grained medical Q&A evaluation**
   Evaluate each item with a prespecified error taxonomy; score omissions and potential harm separately; report results across clinical domain strata; assess confidence calibration; check for benchmark leakage; document adjudication and disagreement-resolution procedures; and version datasets, annotations, rubrics, and benchmark releases for reproducible comparison.

35. **Systematic clinical-note quality benchmark framework**
   Evaluate AI-generated clinical notes across factuality, completeness, hallucination, clinical usefulness, and harm-sensitive scoring; stratify results by note type; report inter-rater reliability; and validate automated metrics against clinician judgments before using them as clinical-quality proxies.

36. **Fine-grained trustworthy domain Q&A evaluation**
   Evaluate domain-specific Q&A datasets with item-level error taxonomies and separate ambiguity and clinically material omission labels; require specialty-expert adjudication; assess calibration and harmful answers; document dataset construction, provenance, annotations, and limitations; and report factual recall separately from clinically safe reasoning.

37. **Clinical-note metric validity and expert-review controls**
   Evaluate factual correctness, completeness, organization, clinical usefulness, hallucination, harmfulness, and style as distinct clinical-note quality dimensions; document each metric's validity and inter-rater reliability; and require expert review when automated scores lack demonstrated validity or agreement with expert judgment for the intended dimension and clinical context.

38. **Blinded intensive-care nursing cross-provider error analysis**
   Benchmark intensive-care nursing question answering across ChatGPT, DeepSeek, and Google Gemini under matched conditions using blinded clinician scoring, and explicitly analyze clinically consequential errors separately from aggregate accuracy.

39. **Trustworthy fine-grained domain QA assessment**
   Evaluate domain-specific medical Q&A with item-level error taxonomies, ambiguity and evidence-quality labels, clinically relevant subgroup analysis, expert adjudication, calibration checks, and dataset documentation covering construction, provenance, annotations, and limitations.

40. **Professional-society guideline concordance evaluation**
   Map each model recommendation to the applicable professional-society statement; score omissions and contradictions separately; and flag advice that exceeds or misstates the guideline's scope.

41. **Plausibility-sensitive clinical answer adjudication**
   Score correctness, clinically important omissions, and risk of harm as separate dimensions; report errors weighted by clinical severity; and require adjudicators to label factually plausible answers as unsafe or incomplete when omitted content or likely downstream use creates clinical risk.

42. **AI-generated clinical-note quality benchmark dimensions**  
   For systematic reviews or experimental benchmarks of AI-generated clinical notes, evaluate factual correctness, omissions, hallucinations, clinical risk, style and format compliance, and inter-rater reliability as separate dimensions; use benchmark datasets with human clinician adjudication to compare evaluation methods.

43. **Clinical-note rubric and safety regression design**  
   For AI-generated clinical-note evaluation, translate systematic review and experimental benchmark findings into prespecified rubrics that separate factuality checks, clinically material omissions, note-quality dimensions, inter-rater adjudication, and safety risk; use the same checks as regression tests before model, prompt, template, retrieval, or integration changes.

44. **Fine-grained specialty-calibrated medical Q&A dataset evaluation**  
   For domain-specific medical Q&A datasets supporting trustworthy medical language models, use an item-level taxonomy; label omissions and potential harms separately; calibrate results by clinical specialty; run adversarial ambiguity and answerability checks; and report benchmark construction, annotation decisions, item-level findings, limitations, and trustworthy-use boundaries.

45. **Domain-specific Q&A benchmark suitability review**  
   Before using a medical Q&A dataset as a trustworthy language-model benchmark, require documented dataset provenance and source-evidence links; specialty-stratified item coverage; an item-level error taxonomy for factual errors, hallucinated or unsupported content, ambiguity, and answerability; separate hallucination and omission scoring; calibration checks by item type or specialty; and explicit suitability criteria that state when the dataset is reliable enough for benchmark reporting.

46. **AI-generated clinical-note evaluation method benchmark reporting**  
   For systematic review-informed experimental benchmarks of AI-generated clinical-note quality methods, compare rubric-based ratings, clinician-rated judgments, automated metrics, factuality checks, omission scoring, and harm-risk evaluations as distinct methods; report benchmark fields for note task and type, dataset provenance, comparator systems, evaluation method, evaluator type, rubric or metric definition, factuality, omission, and harm-risk endpoints, reliability or agreement results, and method-comparison findings.

47. **Fine-grained domain-specific medical Q&A trustworthiness reporting**  
   For trustworthy medical LLM benchmarking, label each domain-specific Q&A item with answerability, ambiguity, evidence support, omission, harm, and uncertainty fields; stratify results by medical subdomain; assess uncertainty calibration within those strata; and document dataset construction, provenance, annotation rules, limitations, and trustworthy-use boundaries.

48. **Ophthalmology CME-style evaluation**  
   For ophthalmology continuing medical education question benchmarks, score correctness, content omission, and risk of harm as separate endpoints; check response consistency with applicable guidelines or reference standards; and report results using specialty-stratified summaries rather than pooled accuracy alone.

49. **ICU nursing scenario question-answering slice**  
   Compare frontier and open models on clinician-authored intensive-care nursing scenarios under matched conditions; report clinically important omissions and unsafe recommendations separately from accuracy; and require human clinical review before any deployment decision.

50. **Trustworthy medical Q&A benchmark construction checks**  
   For domain-specific medical Q&A datasets, construct benchmarks with documented item selection, domain coverage checks, item-level taxonomy labels, and separate omission and hallucination annotations; report annotation practices, dataset limits, and trustworthy-use boundaries alongside model performance.

51. **Stratified domain-specific medical Q&A reporting**  
   For domain-specific medical Q&A evaluation, stratify results by specialty, question type, evidence source, omission risk, calibration, and harm class; compare model performance within those strata and treat aggregate accuracy as a secondary summary rather than the sole trustworthiness claim.

52. **Ophthalmology CME omission-harm separation**  
   For ophthalmology CME-style benchmarks, require scoring schemes that distinguish harmless incompleteness from omissions or recommendations that create clinical risk; report correctness, content omission, and risk of harm as separate dimensions.

53. **Fine-grained domain-specific medical QA grading pattern**  
   For domain-specific medical Q&A evaluation, define question taxonomies, assign separate omission and harm labels, grade answers against supporting evidence, check annotator agreement, analyze clinically relevant subgroups, and report calibration rather than relying on aggregate accuracy alone.

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
- PubMed: Kabir R, Braud SC, Hinson CS, Nazerali RS. "Are large language models consistent with the ASPS and AAPS guidelines? A comparison of AI chatbot recommendations and plastic surgery clinical guidance." *J Plast Reconstr Aesthet Surg*. 2026 May. https://pubmed.ncbi.nlm.nih.gov/41985209/
