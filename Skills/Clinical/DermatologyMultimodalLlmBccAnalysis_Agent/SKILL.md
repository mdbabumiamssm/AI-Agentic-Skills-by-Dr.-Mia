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
name: 'dermatology-multimodal-llm-bcc-analysis'
description: 'Guide multimodal LLM-assisted review of clinical and dermoscopic images for basal cell carcinoma mimickers, emphasizing differential diagnosis, harm checks, and clinician review.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Dermatology Multimodal LLM BCC Analysis

## Overview

Use this skill to support structured multimodal analysis of suspected basal cell carcinoma (BCC) and common mimickers from clinical and dermoscopic images. It is grounded in comparative research on ChatGPT, Gemini, and Claude for BCC image analysis and is designed to improve disciplined differential diagnosis, triage wording, risk checks, and escalation to clinician review.

## When to Use This Skill

- Reviewing clinical photographs or dermoscopic images where BCC is in the differential.
- Comparing BCC with common mimickers such as melanocytic, inflammatory, benign adnexal, or other non-BCC lesions.
- Drafting a clinician-facing differential diagnosis summary from multimodal image observations.
- Producing triage language that avoids definitive diagnosis and emphasizes in-person dermatology assessment.
- Checking whether an AI-generated dermatology image interpretation could create harm through overconfidence, missed melanoma warning signs, or delayed care.
- Preparing structured notes for clinician review, referral, biopsy consideration, or follow-up planning.

## Core Capabilities

1. Image-context intake: Capture lesion site, duration, symptoms, patient risk factors, prior treatments, image modality, and image quality limitations before interpretation.
2. Clinical morphology review: Describe visible features such as pearliness, translucency, ulceration, telangiectasia, pigmentation, scale, border pattern, and anatomic context without overstating certainty.
3. Dermoscopic feature review: Identify reported or visible dermoscopic clues relevant to BCC and mimickers, while flagging low-quality or non-diagnostic images.
4. Differential diagnosis construction: Present BCC and plausible mimickers as ranked or unranked possibilities with supporting and opposing observations.
5. Risk-of-harm screening: Explicitly check for melanoma concern, rapidly changing lesions, bleeding, ulceration, immunosuppression, facial or periocular location, and other features that warrant prompt clinician evaluation.
6. Triage and escalation wording: Generate cautious language for referral urgency, biopsy discussion, follow-up, and patient-facing explanations without claiming a final diagnosis.
7. Clinician review packaging: Produce concise, auditable output that separates observations, uncertainty, differential diagnosis, next-step considerations, and safety caveats.

## Inputs / Outputs

Inputs:

- Clinical lesion image, dermoscopic image, or both.
- Patient-provided or clinician-provided context, including lesion location, duration, symptoms, growth pattern, bleeding, prior skin cancer, immunosuppression, sun exposure history, and relevant medications when available.
- Task intent, such as differential diagnosis, triage support, second-look review, or quality check of an existing AI interpretation.
- Any available clinical constraints, such as teledermatology setting, biopsy access, or need for urgent referral wording.

Outputs:

- Structured image quality and context assessment.
- Morphologic and dermoscopic observations separated from interpretation.
- Differential diagnosis list including BCC and relevant mimickers, with uncertainty clearly stated.
- Risk-of-harm checklist and escalation triggers.
- Clinician-facing summary suitable for review, referral support, or documentation.
- Patient-facing caution statement when requested, emphasizing that AI image review does not replace dermatology evaluation or histopathology when clinically indicated.

## References

- PubMed: ChatGPT, Gemini, and Claude in clinical and dermoscopic image analysis of basal cell carcinoma and its common mimickers: A comparative performance analysis. https://pubmed.ncbi.nlm.nih.gov/41952838/
