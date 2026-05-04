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

- Wang R, Yu C, Dong Q, Qiu J, Wen T. **PKFAR: psychiatry knowledge-fused augmented reasoning with large language models.** PubMed: https://pubmed.ncbi.nlm.nih.gov/41982804/
