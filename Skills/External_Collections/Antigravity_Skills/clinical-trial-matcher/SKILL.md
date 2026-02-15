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

---
name: 'clinical-trial-matcher'
description: 'Matches patient profiles to open clinical trials using vector similarity and inclusion/exclusion criteria. Use when a user provides patient data and asks for relevant trials.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Clinical Trial Matcher

This skill matches a patient's clinical profile against a database of active clinical trials (ClinicalTrials.gov).

## When to use this skill
- When the user provides a patient case summary, medical history, or specific diagnosis.
- When the user explicitly asks to "find trials" or "recruit patients".
- When evaluating eligibility criteria for a specific trial.

## How to use it
1.  **Extract Patient Data:** Identify key demographics (age, sex), primary diagnosis, stage of disease, and prior lines of therapy from the user's input.
2.  **Search Vector Database:**
    -   Generate an embedding for the patient summary.
    -   Query the `trials_vector_db` for top 20 semantic matches.
3.  **Filter by Logic:**
    -   Apply hard filters: `status=Recruiting`, `phase!=Terminated`.
    -   Check inclusion/exclusion criteria (e.g., "Must have EGFR mutation", "No prior immunotherapy").
4.  **Format Output:**
    -   Present matches in a structured list: Trial ID (`NCT number`), Title, Phase, and "Why it fits" (matching logic).
    -   Highlight any potential exclusion criteria that need verification.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->