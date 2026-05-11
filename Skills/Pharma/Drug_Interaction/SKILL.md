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
name: 'drug-interaction-checker'
description: 'Checks for potential drug-drug interactions (DDIs) between a list of medications.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---


# Drug-Drug Interaction (DDI) Checker

This skill analyzes a list of medications to identify known interactions, focusing on safety and contraindications.

## When to Use This Skill

*   Reviewing patient medication lists.
*   Prescribing new medications.
*   Pharmacovigilance monitoring.

## Core Capabilities

1.  **Interaction Detection**: Identifies pairs of drugs with known interactions.
2.  **Severity Grading**: Classifies interactions as Minor, Moderate, or Major.
3.  **Clinical Recommendations**: Provides actionable advice (e.g., "Monitor K+ levels").
4.  **Authoritative Source Comparison Protocol**: When evaluating LLM-generated DDI outputs, compare each interaction, mechanism, recommendation, and severity/category label against authoritative references such as Lexicomp and Drugs.com. For antiseizure medication DDI benchmarking, treat Lexicomp as the primary comparator/reference hierarchy and Drugs.com as an additional comparator, harmonize differing severity/category schemes before comparing outputs, identify high-risk interaction classes from source-defined contraindication/major/severe categories, document antiseizure medication examples separately when present, predefine and report any iterative prompting limit rather than prompting until the model matches a reference, cite the authoritative DDI source databases used rather than relying on model recall, and escalate high-risk combinations for pharmacist or qualified clinician review before using any LLM-only DDI assessment in clinical decision support.
5.  **LLM DDI Evidence Hierarchy and Validation Checklist**: For LLM-generated DDI answers, rank curated references such as Lexicomp or Drugs.com above model-generated content, use antiseizure medication DDIs as a stress case when applicable, validate whether the answer identifies the interacting drug pair, mechanism, severity/category, clinical recommendation, and source agreement/disagreement, and retain iterative prompt audit logs showing each prompt, answer revision, and comparison outcome against the curated reference.

## Workflow

1.  **Input**: List of drug names (e.g., "Warfarin, Aspirin").
2.  **Analysis**: Queries internal interaction database.
3.  **Output**: Interaction report with severity and mechanisms.

## Example Usage

**User**: "Check interactions for Warfarin and Aspirin."

**Agent Action**:
```bash
python3 Skills/Pharma/Drug_Interaction/impl.py --drugs "Warfarin, Aspirin"
```

## References

*   https://pubmed.ncbi.nlm.nih.gov/41994367/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
