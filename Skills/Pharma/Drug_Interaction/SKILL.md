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
10. **LLM-Based DDI Benchmark Caution**: For LLM-based DDI checks, use antiseizure medication interactions as a high-risk benchmark example; require source grounding against authoritative references such as Lexicomp or product labeling, reproducible iterative prompt checks, severity classification, and pharmacist review before clinical use.
11. **Lexicomp-Centered LLM DDI Benchmarking**: When benchmarking LLM drug-drug interaction responses or consumer reference outputs for antiseizure medications, compare each response against an authoritative reference such as Lexicomp, set a finite iterative prompting limit before testing, track whether each iteration changes the interaction pair, severity/category, mechanism, recommendation, or source agreement, and escalate to pharmacist or qualified clinician review whenever Lexicomp, Drugs.com-style references, and LLM outputs disagree.
12. **ASM DDI Reference Comparison and Citation Control**: For antiseizure medication DDI evaluation, compare LLM responses with curated references such as Lexicomp and Drugs.com, normalize interaction severity labels before reporting agreement or disagreement, cite the reference used for each clinically relevant interaction, state the bounded iterative prompting limit, and require pharmacist oversight before clinical use.
13. **Lexicomp/Drugs.com-Anchored ASM Benchmark Workflow**: For antiseizure medication DDI benchmarking, anchor comparisons to Lexicomp and Drugs.com where available, predefine and report the iterative prompting limit, check agreement in severity-class assignments before accepting outputs, require citations for each clinically relevant interaction, and keep pharmacist review mandatory before any clinical use.
14. **LLM DDI Confidence Threshold Gate**: For LLM-based drug-drug interaction checks, benchmark answers against trusted references such as Lexicomp and Drugs.com, treat antiseizure medication interactions as edge cases requiring extra scrutiny, predefine confidence thresholds for accepting, flagging, or rejecting model outputs, treat iterative prompting as a potential answer-drift risk rather than proof of correctness, and require pharmacist review before any clinical recommendation is used.
15. **Antiseizure DDI Comparative Review**: For antiseizure medication DDI evaluations, treat Lexicomp as the primary comparator and Drugs.com as a secondary comparator when available; normalize severity labels before judging agreement with LLM output, document any iterative prompting steps as exploratory caveats rather than validation, and require pharmacist review before relying on the final interaction assessment.
16. **ASM Interaction Prompt-Sensitivity Safeguards**: For antiseizure medication interaction checks, use an authoritative compendium as the benchmark reference, test and report whether iterative prompting changes the result, evaluate interaction severity separately from management advice, cite the supporting evidence, and never offer reassurance produced through iterative prompting unless the benchmark evidence supports it.
17. **Source-Dated ASM DDI Benchmarking**: Benchmark antiseizure-medication DDI outputs against authoritative references such as Lexicomp, recording each source and access/update date; use a predefined, finite iterative-prompting protocol; extract interaction severity and mechanism separately; document disagreements among the LLM, Drugs.com, and Lexicomp without resolving them by model consensus; and escalate unresolved, clinically significant, or source-date-sensitive findings to a pharmacist before clinical use.
18. **Validated ASM DDI Evaluation**: Evaluate antiseizure-medication DDI outputs against authoritative references such as Lexicomp, test sensitivity to iterative prompting, verify interaction severity and mechanism, check for omitted interactions, and never use LLM output as a replacement for a validated interaction database.

## Benchmark Caution: LLM-Based DDI Checks

LLM-generated DDI findings are screening drafts, not standalone clinical evidence. For antiseizure medication interaction checks, ground each claimed interaction in authoritative references such as Lexicomp or product labeling, document prompts and iterations so another reviewer can reproduce the result, classify severity explicitly, and route outputs to pharmacist review before use in patient care.

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
