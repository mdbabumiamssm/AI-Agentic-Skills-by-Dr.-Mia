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
4.  **LLM-vs-Reference DDI Checking Caution**: For antiseizure medication interaction checks, use LLM output only as a draft that must be verified against Lexicomp/Drugs.com-style DDI references. Reconcile severity/category differences across sources, document any finite iterative prompting limit and remaining disagreements, and route unresolved or clinically significant interactions to pharmacist or qualified clinician review before clinical use.
5.  **LLM DDI Evidence Hierarchy and Validation Checklist**: For LLM-generated DDI answers, rank curated references such as Lexicomp or Drugs.com above model-generated content, use antiseizure medication DDIs as a stress case when applicable, validate whether the answer identifies the interacting drug pair, mechanism, severity/category, clinical recommendation, and source agreement/disagreement, and retain iterative prompt audit logs showing each prompt, answer revision, and comparison outcome against the curated reference.
6.  **Antiseizure Medication DDI Evaluation**: For antiseizure medication interaction checks, compare LLM outputs with Lexicomp and Drugs.com references, reconcile interaction severity/category differences before presenting conclusions, set and report a finite iterative prompting limit, and require pharmacist or qualified clinician verification before clinical use.
7.  **Authoritative ASM DDI Review**: When evaluating antiseizure medication DDIs, classify interaction severity from authoritative references such as Lexicomp, cite the evidence source used for each clinically relevant interaction, apply iterative prompting only as a bounded safeguard to improve LLM output, and require pharmacist review before recommendations are used for patient care.
8.  **ASM DDI Benchmarking Guidance**: When benchmarking antiseizure medication DDI checks, compare outputs against authoritative references such as Lexicomp and Drugs.com, audit severity agreement and missing interactions, treat iterative prompting as a risk surface that may change answers without proving correctness, and require pharmacist review for unresolved, missing, or clinically significant interaction findings.
9.  **ASM DDI Benchmark Source Priority**: For antiseizure medication DDI benchmarking against Lexicomp and Drugs.com, treat curated references as higher priority than LLM outputs, define the iterative prompting cap before review, reconcile severity discrepancies explicitly, and warn users whenever an LLM answer conflicts with curated reference findings.

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
