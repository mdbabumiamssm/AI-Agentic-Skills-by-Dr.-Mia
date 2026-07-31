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
name: 'psychiatry-knowledge-fused-reasoning'
description: 'Guides psychiatry-specific LLM support using PKFAR-style knowledge fusion and augmented reasoning for diagnosis, risk, medications, and care planning.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Psychiatry Knowledge-Fused Reasoning

## Overview

This skill guides psychiatry-focused LLM decision support using a knowledge-fused augmented reasoning pattern inspired by PKFAR. It helps structure psychiatric reasoning around evidence-grounded clinical facts, longitudinal symptoms, safety risks, differential diagnosis, medications, and guideline-sensitive recommendations while preserving clinical uncertainty.

The skill is intended for clinician-facing support, chart review, education, and structured reasoning assistance. It must not be used as a standalone diagnostic or treatment authority.

## When to Use This Skill

- The task involves psychiatric diagnostic reasoning, differential diagnosis, or symptom formulation.
- The user asks for psychiatry-specific clinical decision support grounded in domain knowledge.
- The case includes longitudinal psychiatric symptoms, functional impairment, comorbidity, medication history, or psychosocial context.
- The work requires explicit handling of suicide risk, violence risk, self-neglect, substance use, psychosis, mania, trauma, or medication safety.
- The output should separate observed facts, patient-reported history, inferred hypotheses, evidence gaps, and recommended next clinical questions.
- The user needs a structured psychiatric reasoning artifact for review, teaching, triage, or care-planning support.

## Core Capabilities

1. **Knowledge-grounded case framing**  
   Organize psychiatric information into clinically meaningful domains: presenting concern, symptom timeline, mental status findings, risk factors, protective factors, medications, substance use, medical contributors, psychosocial context, and prior treatment response.

2. **Augmented reasoning trace**  
   Build a transparent reasoning path that distinguishes documented facts from hypotheses, flags uncertainty, and identifies missing data needed before strong conclusions.

3. **Differential diagnosis support**  
   Generate and compare plausible psychiatric and medical differentials using symptom clusters, time course, impairment, exclusions, and comorbidities without prematurely collapsing uncertainty.

4. **Risk and safety synthesis**  
   Summarize acute and chronic risk considerations for self-harm, suicide, violence, grave disability, withdrawal, intoxication, medication adverse effects, and need for higher level of care.

5. **Medication and treatment reasoning**  
   Review medication history, adherence, response, adverse effects, interactions, contraindications, and monitoring needs while keeping recommendations guideline-sensitive and clinician-verifiable.

6. **Clinical communication**  
   Produce structured outputs for clinicians, including concise assessment, problem representation, next questions, evidence gaps, and suggested care-planning considerations.

7. **PKFAR-style knowledge fusion checks**  
   Retrieve curated psychiatric knowledge bases when available, use augmented reasoning steps around case evidence, run safety checks for diagnosis, medication, and risk assessment, report uncertainty explicitly, and limit outputs to clinician-overseen support rather than standalone psychiatric authority.

8. **PKFAR knowledge-fusion workflow and evaluation**  
   Apply explicit retrieval and fusion stages: retrieve psychiatry-relevant evidence, attribute each fused claim to its supporting case fact or source, and apply psychiatry-specific reasoning controls before synthesis. Preserve competing differential diagnoses with supporting, opposing, and missing evidence; require medication safeguards for contraindications, interactions, adverse effects, monitoring, and clinician verification; and require structured suicide, self-harm, violence, intoxication, withdrawal, and level-of-care risk checks with urgent escalation when indicated. When evaluating an implementation, use ablation testing to assess the contribution of retrieval, fusion, and reasoning controls, and use clinician-adjudicated review rather than unsupported automated correctness claims.

9. **PKFAR core psychiatry pattern**  
   Use PKFAR as the core psychiatry knowledge-fused reasoning pattern: retrieve only case-relevant psychiatric knowledge from user-provided records, local policies, guidelines, formularies, medication references, and trusted clinical sources; fuse retrieved knowledge with the supplied case context; triage suicide, self-harm, violence, psychosis, mania, intoxication, withdrawal, self-neglect, and grave-disability risks; check medication decisions against diagnosis, comorbidity, substance use, pregnancy, labs, adverse effects, interactions, adherence, prior response, and monitoring context; state uncertainty, missing data, and competing explanations; and escalate urgent or high-risk findings to clinician-supervised emergency or higher-level-of-care review instead of issuing autonomous decisions.

10. **PKFAR motivating architecture**  
   Use PKFAR as the motivating architecture for psychiatry support: combine knowledge-fused retrieval with reasoning augmentation, run diagnostic and medication safety checks, enforce risk assessment constraints, and document which statements come from case evidence or retrieved sources versus model inference.

11. **PKFAR citation and escalation guardrails**  
   Use structured psychiatric knowledge retrieval before synthesis; cite the supplied case fact, retrieved psychiatric source, or model-inference status for diagnostic, medication, and risk-relevant conclusions; avoid final diagnosis, prescribing, dosing, or disposition decisions; and escalate urgent, high-risk, or clinically uncertain findings to a licensed clinician or emergency care pathway.

12. **PKFAR-style clinician-reviewed support loop**  
   Use structured retrieval from psychiatry knowledge sources to support symptom and medication reasoning, risk assessment, and provenance-labeled synthesis; route outputs through clinician review for mental health decision support, especially when findings affect diagnosis, medication changes, safety planning, or level-of-care decisions.

13. **PKFAR psychiatry-specific knowledge injection pattern**  
   Apply PKFAR as a psychiatry-specific knowledge-fused augmented reasoning pattern by retrieving and injecting relevant psychiatric knowledge before synthesis; enforce diagnostic-risk controls, medication safety checks, hallucination checks against case facts and retrieved sources, and clinician-reviewed care planning rather than autonomous diagnosis, prescribing, or disposition.

14. **PKFAR-style evidence-bound psychiatry synthesis**  
   Use structured psychiatric knowledge retrieval before reasoning; keep diagnosis and risk conclusions within the limits of supplied case facts and retrieved sources; ground medication and care-plan considerations in cited evidence, contraindications, monitoring needs, and clinician review; display uncertainty, missing data, and competing explanations; and block unsupported psychiatric recommendations, including autonomous diagnosis, prescribing, dosing, disposition, or safety-plan decisions.

15. **PKFAR-style structured retrieval and oversight pattern**  
   Apply a structured retrieval pass before psychiatric synthesis; inject diagnostic, risk, and medication knowledge only when it is case-relevant and attributable to supplied facts or retrieved sources; label evidence provenance for clinically important claims; and keep final diagnosis, prescribing, dosing, disposition, and safety-plan decisions under human psychiatric oversight.

16. **PKFAR psychiatry-specific source-vetted reasoning pattern**  
   Apply PKFAR as a psychiatry-specific knowledge-fused reasoning pattern by vetting retrieval sources before use, supporting psychiatric differential diagnosis without replacing clinical judgment, checking medication and safety-risk implications, structuring uncertainty around missing or conflicting evidence, and requiring clinician review for any diagnosis, medication, safety, or disposition-relevant output.

## Inputs / Outputs

**Inputs**

- De-identified psychiatric history, clinical note, transcript, intake form, or structured case summary.
- Presenting symptoms, timeline, prior diagnoses, current and past medications, treatment response, substance use, medical history, family history, trauma history, and social context.
- Mental status exam findings, risk assessment data, collateral information, laboratory or medical findings, and level-of-care constraints when available.
- Local guideline, formulary, policy, or institution-specific reference material when the user provides it.

**Outputs**

- A structured psychiatric formulation grounded in the supplied case facts.
- A differential diagnosis table or ranked list with supporting and opposing evidence.
- A risk synthesis that separates acute risk, chronic risk, protective factors, warning signs, and missing safety data.
- Medication and treatment considerations with cautions, monitoring needs, and clear clinician-verification requirements.
- Follow-up questions and recommended information to obtain before making or changing a clinical plan.
- A brief safety note when the case suggests imminent risk, emergency evaluation, mandated reporting, or urgent clinician escalation may be needed.

## References

- Wang R, Yu C, Dong Q, Qiu J, Wen T. **PKFAR: psychiatry knowledge-fused augmented reasoning with large language models.** Health Inf Sci Syst, 2026 Dec. PubMed PMID: 41982804. https://pubmed.ncbi.nlm.nih.gov/41982804/
