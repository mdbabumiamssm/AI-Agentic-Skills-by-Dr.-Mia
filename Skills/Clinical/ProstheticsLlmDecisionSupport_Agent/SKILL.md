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
name: 'prosthetics-llm-decision-support'
description: 'Assist clinician-reviewed prosthetic recommendations using ProsthetiX-AI-style intake, evidence retrieval, safety checks, and documentation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Prosthetics LLM Decision Support

## Overview

Use this skill to support evidence-based prosthetic recommendation workflows for people with limb loss, especially when an agent must organize intake, retrieve supporting evidence, screen for safety constraints, and prepare clinician-reviewable options. The skill is assistive only: it must not prescribe, authorize coverage, replace prosthetist judgment, or present recommendations as autonomous clinical decisions.

## When to Use This Skill

- A user asks for prosthetic component, socket, suspension, liner, interface, or rehabilitation recommendation support.
- A case includes amputation level, functional status, comorbidities, residual-limb issues, patient goals, or payer constraints that need structured decision support.
- A clinician-facing note, options table, evidence summary, or prior-authorization rationale is needed for prosthetic planning.
- A workflow needs contraindication, compatibility, weight-limit, skin-integrity, infection-risk, fall-risk, or material-allergy checks before prosthetic options are discussed.
- A patient-facing explanation is needed after clinician review, using plain language and clear uncertainty.
- Do not use this skill as the sole basis for urgent wound, infection, ischemia, severe pain, post-operative complication, or emergency-care decisions.

## Core Capabilities

1. Intake normalization: Collect and structure amputation level, laterality, residual-limb shape and skin status, time since amputation, height and weight, mobility classification, gait aids, activity level, home and work environment, cognition, dexterity, prior prosthesis history, pain, goals, and clinician constraints.

2. Functional and goal alignment: Map candidate options to confirmed or clinician-estimated functional level, rehabilitation goals, terrain needs, occupational demands, fall history, endurance, and expected training capacity without assigning a definitive functional level from narrative alone.

3. Evidence retrieval: Use current PubMed, clinical guidance, manufacturer documentation, payer policy, and local formulary sources when available. Prefer primary or authoritative sources, date-stamp retrieved evidence, and separate evidence-supported statements from inference.

4. Safety and contraindication screening: Check wounds, edema, infection concern, diabetes or vascular disease, neuropathy, skin fragility, high fall risk, weight or activity limits, component compatibility, material sensitivity, waterproofing needs, battery or charging constraints, and cognitive or manual ability to operate advanced components.

5. Recommendation framing: Present options as clinician-reviewable alternatives with rationale, uncertainties, missing data, expected benefits, tradeoffs, and follow-up needs. Include conservative alternatives when high-cost or advanced components are considered.

6. Payer and documentation support: Identify documentation elements commonly needed for coverage review, such as functional classification, medical necessity, weight limit, trial history, expected use environment, and therapy plan. Do not guarantee reimbursement or encourage upcoding.

7. Clinician handoff: End with a concise review checklist, unresolved questions, evidence links, and explicit language that the output requires prosthetist or prescribing clinician review before clinical use.

8. ProsthetiX-AI example workflow: For evidence-based prosthetic recommendations, sequence structured prosthetic intake, device constraint capture, evidence retrieval, functional goal mapping, recommendation rationale, safety and contraindication checks, cited recommendation summaries, shared decision documentation, and explicit clinician approval before any patient-facing advice is finalized or shared.

9. ProsthetiX-AI workflow details: Build a patient intake schema covering de-identified demographics, amputation details, residual-limb condition, comorbidities, current prosthesis, mobility/function goals, activity environment, preferences, and device constraints; run contraindication checks and safety review before option ranking; retrieve evidence before drafting recommendation rationale; and output structured intake, evidence summary, safety flags, recommendation rationale, unresolved questions, and documentation for prosthetist oversight and sign-off.

## Inputs / Outputs

Inputs this skill can consume:

- Patient case summary with de-identified demographics, amputation level, limb condition, comorbidities, functional status, goals, and current prosthesis history.
- Clinician findings, therapy notes, gait assessment, residual-limb measurements, skin observations, and weight-bearing tolerance.
- Payer, region, plan type, local formulary, manufacturer, or component constraints if provided.
- User preference for clinician-facing, patient-facing, or payer-documentation output.

Outputs this skill should produce:

- Structured intake summary with missing critical data highlighted.
- Safety and contraindication checklist, including urgent issues that need clinical escalation.
- Evidence table with citation links, source type, retrieval date when web evidence is used, and relevance to the case.
- Ranked or grouped prosthetic options covering socket/interface, suspension, liner, foot/ankle, knee, terminal device, or adjunctive rehabilitation needs as applicable.
- Rationale that distinguishes evidence, policy constraints, manufacturer limits, clinician judgment points, and LLM inference.
- Clinician-review note, payer documentation draft, or patient-facing explanation that states the output is not a final prescription or coverage determination.

## References

- Kumar V, Pratihar DK. "ProsthetiX-AI: An LLM-based clinical decision support system for evidence-based prosthetic recommendations." PubMed PMID: 41978836. https://pubmed.ncbi.nlm.nih.gov/41978836/
