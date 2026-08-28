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
name: 'surgical-patient-communication-ai'
description: 'Draft clinician-reviewed surgical communication, education, consent-support, expectation-setting, and escalation materials using AI governance for surgery workflows.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Surgical Patient Communication AI

## Overview

Use this skill to create patient-facing and clinician-facing communication artifacts for surgical care, including preoperative education, informed-consent support, perioperative expectation setting, postoperative instructions, and escalation guidance. The workflow is grounded in the finding that artificial intelligence tools may improve patient-physician communication in surgery while requiring careful governance, clinician review, and clear boundaries.

This skill supports communication drafting and review; it does not replace the surgeon, consent process, diagnosis, treatment planning, or institution-specific legal and clinical requirements.

## When to Use This Skill

- Draft plain-language explanations of a surgical procedure, indication, expected course, alternatives, benefits, risks, and uncertainties.
- Prepare clinician-reviewed informed-consent support materials without presenting the AI output as consent itself.
- Convert surgeon-approved perioperative instructions into patient-friendly summaries, FAQs, scripts, checklists, or portal-message drafts.
- Adapt surgical education materials for health literacy, language access, cultural sensitivity, accessibility, or caregiver involvement.
- Build teach-back prompts, patient question lists, shared decision-making aids, or expectation-setting materials for surgery.
- Identify communication gaps, missing safety warnings, red-flag symptoms, or escalation points that require a clinician.
- Review AI-generated surgical communication for unsafe reassurance, missing risks, scope creep, unsupported claims, or lack of clinician handoff.

## Core Capabilities

1. **Procedure-specific communication drafting**  
   Create concise, plain-language explanations of what will happen before, during, and after surgery using the procedure details supplied by the clinical team.

2. **Informed-consent support governance**  
   Help organize risks, benefits, alternatives, expected recovery, and uncertainty while preserving the requirement that a qualified clinician conduct and document consent.

3. **Perioperative education and expectation setting**  
   Produce preoperative preparation instructions, day-of-surgery guidance, recovery timelines, activity restrictions, medication reminders, wound-care explanations, and follow-up expectations from approved source material.

4. **Teach-back and comprehension checking**  
   Generate patient-centered questions and prompts that help clinicians confirm understanding, surface concerns, and correct misconceptions.

5. **Escalation and safety boundary detection**  
   Flag symptoms, decisions, emotional distress, urgent postoperative concerns, or procedure-specific complications that should be routed to a clinician or emergency pathway.

6. **Health-literacy and accessibility adaptation**  
   Rewrite materials in plain language, avoid unexplained jargon, preserve essential clinical meaning, and include accommodations for caregivers, interpreters, sensory needs, or limited literacy when requested.

7. **Clinician review package creation**  
   Return drafts with assumptions, missing information, safety flags, and suggested clinician-review checkpoints so the output can be verified before patient use.

## Inputs / Outputs

### Inputs

- Surgical procedure name, indication, laterality or site, urgency, planned setting, and perioperative phase.
- Surgeon-approved or institution-approved patient education, consent language, order sets, discharge instructions, or clinic templates.
- Patient context relevant to communication, such as age group, preferred language, literacy needs, caregiver role, disability accommodation, anxiety concerns, or common misconceptions.
- Known clinical details that the treating team has authorized for communication, including diagnosis, relevant comorbidities, medication restrictions, anesthesia plan, and recovery constraints.
- Local escalation pathways, emergency instructions, contact numbers, after-hours workflow, and clinician review requirements.
- Requested output type, such as script, handout, portal reply, FAQ, checklist, visit agenda, teach-back guide, or clinician handoff note.

### Outputs

- Patient-facing draft communication that is clear, respectful, actionable, and bounded by the supplied clinical facts.
- Clinician-facing review notes listing assumptions, missing details, safety concerns, and places where local policy or surgeon judgment is required.
- Informed-consent support outline covering procedure purpose, expected benefits, material risks, alternatives, no-treatment option, recovery expectations, and patient questions.
- Teach-back prompts and comprehension checks tailored to the surgical scenario.
- Escalation guidance that separates routine questions from urgent symptoms and clinician-only decisions.
- A final review reminder that the output requires clinician approval before use in patient care.

## References

- Clement EA, Lee L. Artificial Intelligence to Improve Patient-Physician Communication in Surgery. *Clinics in Colon and Rectal Surgery*. 2026 May. PMID: 41948164. https://pubmed.ncbi.nlm.nih.gov/41948164/
- Zhou XY, Guo Y, Shen M, Yang GZ. Artificial Intelligence in Surgery. arXiv:2001.00627. https://arxiv.org/abs/2001.00627
- Tu T, Palepu A, Schaekermann M, et al. Towards Conversational Diagnostic AI. arXiv:2401.05654. https://arxiv.org/abs/2401.05654
