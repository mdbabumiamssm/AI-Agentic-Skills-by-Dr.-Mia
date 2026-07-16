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
name: 'gastroenterology-llm-clinical-practice'
description: 'Support evidence-grounded use of LLMs in gastroenterology clinical practice, including endoscopy, IBD, hepatology, triage, documentation, and patient education.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Gastroenterology LLM Clinical Practice

## Overview

Use this skill to structure safe, evidence-grounded large language model support for gastroenterology clinical practice workflows. It helps route requests involving endoscopy explanation, inflammatory bowel disease, hepatology, documentation, triage, and patient education while preserving clinician oversight and avoiding autonomous diagnosis or treatment.

This skill is grounded in a 2026 scoping review of LLM applications in gastroenterology clinical practice and should be used to make outputs clinically cautious, source-aware, and scoped to decision support.

## When to Use This Skill

- Use when a user asks for LLM-supported gastroenterology clinical workflow design, implementation, or evaluation.
- Use for gastroenterology patient education drafts, discharge instructions, procedure explanations, bowel preparation guidance, or plain-language summaries.
- Use for clinician-facing documentation support, such as note drafting, referral summaries, handoff summaries, or structured chart review in GI contexts.
- Use for endoscopy-adjacent support, including patient-facing procedure explanation, findings summarization, or follow-up communication, without interpreting images independently.
- Use for IBD, hepatology, functional GI, motility, pancreaticobiliary, or general GI decision-support prompts when clinical uncertainty must be preserved.
- Use for triage language that helps identify urgency signals and recommends clinician or emergency evaluation when appropriate.
- Use when reviewing or designing safeguards for GI LLM use, including evidence grounding, disclaimers, escalation criteria, privacy checks, and human review.
- Do not use as a substitute for licensed clinical judgment, direct diagnosis, medication prescribing, procedure selection, or emergency care.

## Core Capabilities

1. Scope the clinical task.
   Identify whether the request is patient education, documentation, triage support, literature-grounded explanation, workflow design, or quality/safety review. State the intended user, clinical setting, and whether the output is patient-facing or clinician-facing.

2. Enforce safety boundaries.
   Avoid autonomous diagnosis, treatment selection, dosing, procedural decisions, or interpretation of raw endoscopy images. Flag red-flag symptoms, unstable presentations, suspected GI bleeding, severe abdominal pain, acute liver failure signs, dehydration, sepsis concern, or other urgent features for immediate clinician or emergency evaluation.

3. Ground outputs in evidence.
   Prefer current clinical guidelines, institutional policies, drug labels, and peer-reviewed sources when available. Distinguish sourced facts from general clinical framing, and state when recommendations require local guideline or clinician confirmation.

4. Support patient communication.
   Convert GI concepts into plain language while preserving clinically important cautions. Use readable structure for bowel preparation, procedure expectations, medication questions, follow-up steps, warning signs, and when to seek urgent help.

5. Support clinician documentation.
   Draft concise, review-ready clinical text such as histories, assessment summaries, referral letters, endoscopy follow-up messages, and patient instructions. Preserve uncertainty and avoid adding facts not present in the source note.

6. Assist GI triage workflows.
   Organize symptoms, duration, severity, risk factors, medication exposures, prior diagnoses, and alarm features. Produce routing suggestions as decision support only, with explicit clinician review requirements.

7. Review LLM implementation risks.
   Check for hallucinated citations, outdated guidance, missing escalation criteria, privacy exposure, unverified patient-specific recommendations, biased assumptions, and absence of human-in-the-loop review.

8. Adapt to gastroenterology subdomains.
   Tailor framing for endoscopy, IBD, hepatology, colorectal screening, functional GI disorders, pancreaticobiliary disease, nutrition, and medication counseling while keeping outputs within the available evidence and requested scope.

## Inputs / Outputs

### Inputs

- Clinical context: patient-facing, clinician-facing, operational, research, or implementation review.
- GI subdomain: endoscopy, IBD, hepatology, colorectal screening, functional GI, pancreaticobiliary, motility, nutrition, or general gastroenterology.
- Source material: notes, referral text, guideline excerpts, patient question, draft prompt, workflow description, or literature summary.
- Intended output type: education handout, triage script, documentation draft, prompt template, safety checklist, workflow review, or evidence-grounded explanation.
- Constraints: reading level, language, local policy, urgent-care thresholds, required citations, privacy requirements, and clinician review process.

### Outputs

- A scoped gastroenterology LLM response or workflow artifact aligned with the requested use case.
- Clear separation between patient-facing language and clinician-facing decision support.
- Evidence-grounding notes or citation placeholders when source verification is required.
- Explicit limitations, uncertainty statements, and clinician oversight requirements.
- Red-flag escalation language when symptoms or findings could require urgent evaluation.
- A safety review checklist for clinical deployment or prompt evaluation when requested.

## References

- PubMed: Yazarkan Y, Sonmez G, Simsek C. "Artificial intelligence in gastroenterology clinical practice: Scoping review of large language model applications." Int J Med Inform. Published 2026 Jul 1. https://pubmed.ncbi.nlm.nih.gov/41921370/
