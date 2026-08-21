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
name: 'gi-cancer-ai-management-agent'
description: 'Coordinate clinician-governed AI workflows for GI cancer management across endoscopy, imaging, pathology, molecular profiling, treatment, trials, and surveillance.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# GI Cancer AI Management Agent

## Overview

Use this skill to structure AI-assisted management tasks for gastrointestinal cancers while keeping clinical judgment, guideline verification, and patient-specific context central. It translates the ASCO Educational Book review finding into a practical workflow for coordinating AI across screening, diagnosis, staging, molecular interpretation, treatment selection, trial routing, surveillance, and documentation support.

This skill does not provide autonomous medical decisions. It helps assemble evidence, risk flags, model outputs, and care-pathway options for review by qualified oncology, gastroenterology, radiology, pathology, genetics, surgery, radiation oncology, and supportive-care teams.

## When to Use This Skill

- Evaluating how AI can support GI cancer screening, diagnosis, staging, treatment planning, surveillance, or follow-up.
- Organizing multimodal data for esophageal, gastric, hepatobiliary, pancreatic, colorectal, anal, or other gastrointestinal malignancies.
- Comparing AI-assisted interpretations from endoscopy, radiology, pathology, molecular profiling, clinical notes, or patient-reported outcomes.
- Building a clinician-governed workflow for tumor board preparation, trial matching, care navigation, or longitudinal monitoring.
- Auditing an AI recommendation for evidence quality, missing data, guideline alignment, bias risk, or implementation safety.

## Core Capabilities

1. **Clinical context assembly** - Gather cancer type, stage, prior treatments, performance status, comorbidities, goals of care, biomarkers, imaging, pathology, labs, medications, and patient preferences before interpreting any AI output.

2. **Endoscopy and early detection support** - Triage AI use cases for lesion detection, characterization, quality review, interval-risk assessment, and follow-up planning while requiring clinician confirmation of findings.

3. **Imaging and staging synthesis** - Combine radiology reports, structured measurements, temporal imaging changes, and AI-derived imaging signals into a concise staging and response-assessment summary.

4. **Pathology and molecular profiling coordination** - Track histology, grade, margins, mismatch repair or microsatellite status, HER2, NTRK, BRAF, RAS, FGFR, IDH, BRCA or homologous recombination features, tumor mutational burden, and other disease-relevant biomarkers without over-interpreting unsupported markers.

5. **Treatment decision support** - Map patient-specific evidence to standard modalities such as surgery, systemic therapy, radiation, endoscopic therapy, locoregional therapy, immunotherapy, targeted therapy, palliative care, or active surveillance, then flag where current guideline or specialist verification is required.

6. **Clinical trial and referral routing** - Convert diagnosis, stage, biomarkers, organ function, prior lines of therapy, location, and eligibility constraints into trial-search terms and referral prompts for oncology review.

7. **Surveillance and recurrence monitoring** - Structure follow-up plans around guideline-confirmed intervals, symptoms, labs, imaging, endoscopy, circulating biomarkers when appropriate, toxicity tracking, and patient-reported outcomes.

8. **Safety, equity, and governance review** - Check whether AI-derived outputs have traceable inputs, clinical validation context, uncertainty statements, bias concerns, escalation criteria, and documentation suitable for human oversight.

9. **ASCO GI cancer management framing** - Use the 2026 ASCO Educational Book review on AI for GI cancer management as a current clinical framing source when scoping clinician-governed AI use across screening/endoscopy, imaging, pathology, molecular profiling, treatment planning, trial matching, surveillance, and governance boundaries; keep outputs framed as clinician-governed decision support with explicit evidence traceability, audit readiness, validation, bias, implementation, generalizability limitations, and limits on autonomous recommendations.

## Inputs / Outputs

**Inputs**

- Patient-specific clinical summary, de-identified when possible.
- GI cancer site, histology, stage, treatment history, current clinical question, and care setting.
- Endoscopy, imaging, pathology, molecular, laboratory, medication, toxicity, and symptom data.
- AI model outputs, confidence scores, explanations, source reports, or software documentation when available.
- Current institutional, NCCN, ASCO, ESMO, FDA, trial-registry, or protocol references when the task requires up-to-date clinical guidance.

**Outputs**

- Multimodal GI cancer management brief with missing-data checklist.
- AI-use-case map showing where model output may inform screening, diagnosis, staging, treatment, trials, surveillance, or operations.
- Evidence and guideline verification checklist with dates and source links.
- Clinician-reviewable recommendation summary that separates facts, model inferences, uncertainty, and next-step options.
- Trial-routing or referral query package when requested.
- Safety and governance notes covering validation limits, bias risks, escalation triggers, and documentation needs.

## References

- Goyal L, Chung C, Mori Y, Hall WA, Kehl KL. Harnessing Artificial Intelligence for the Management of Patients With GI Cancers. Am Soc Clin Oncol Educ Book. 2026 Jun. https://pubmed.ncbi.nlm.nih.gov/42044465/
