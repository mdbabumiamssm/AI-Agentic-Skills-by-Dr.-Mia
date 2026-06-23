---
name: 'drug-interaction-checker'
description: 'Checks for potential drug-drug interactions (DDIs) between a list of medications.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

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
19. **ASM DDI Validation and Review Gates**: Validate antiseizure-medication interaction cases against authoritative drug-information references such as Lexicomp, use predefined iterative-prompting controls, measure and report agreement or disagreement for severity and mechanism separately, and require both a supporting citation and pharmacist review before accepting any clinically relevant output.
20. **ASM DDI Benchmark Procedure**: Benchmark antiseizure-medication interactions against an authoritative reference such as Lexicomp; score severity and mechanism agreement separately; repeat the evaluation with bounded iterative prompts to identify answer sensitivity; verify that every interaction claim has a supporting citation; and mandatorily escalate uncertain, conflicting, or unsupported interactions to a pharmacist or qualified clinician.
21. **Comparative DDI Benchmark Audit**: Compare LLM antiseizure-medication interaction assessments with authoritative compendia such as Lexicomp and Drugs.com; score severity and mechanism agreement separately; test and record sensitivity to bounded iterative prompting; track omitted interactions; verify supporting citations; and escalate conflicting source findings to a pharmacist or qualified clinician.
22. **PubMed 41994367 ASM DDI Benchmarking**: For antiseizure medication DDI benchmarks informed by the 2026 cross-sectional comparison of LLMs and Drugs.com versus Lexicomp, compare generated answers against authoritative databases rather than model consensus, report severity-classification agreement or disagreement, state bounded iterative prompting limits, and require pharmacist review before clinical interpretation.
23. **ASM DDI Omission and Recommendation Audit**: When checking antiseizure medication DDIs, compare LLM outputs against authoritative references such as Lexicomp, capture omissions in severity, mechanism, or evidence support, preserve iterative prompt audit trails, and avoid clinical recommendations that are not supported by the benchmark reference.
24. **Authoritative DDI Evaluation and Prompt Testing**: Prioritize validated interaction databases such as Lexicomp over Drugs.com and LLM output; assess severity and mechanism agreement separately, treat omitted interactions as false-negative risks, and use bounded iterative prompting to test answer stability rather than establish correctness.
25. **Authoritative Benchmark and Escalation Record**: Benchmark antiseizure-medication examples against an authoritative reference such as Lexicomp; run bounded iterative-prompt sensitivity testing; report interaction severity and mechanism as separate fields; and explicitly escalate to a pharmacist or qualified clinician whenever evidence conflicts across sources or is incomplete.
26. **Benchmark-Driven DDI Verification**: Normalize medication names and formulations before comparison; verify each DDI against an authoritative reference such as Lexicomp; report interaction-severity concordance; cite the evidence supporting each clinically relevant finding; predefine and log finite prompt iterations; escalate conflicting, incomplete, or unsupported findings to a pharmacist or qualified clinician; and never treat repeated prompting or repeated model agreement as independent confirmation.
27. **Authoritative ASM Interaction Benchmarking**: Benchmark antiseizure-medication interaction cases against Lexicomp-style authoritative references; test and record prompt sensitivity with bounded iterative prompting; report severity concordance; verify every supporting citation; and mandatorily escalate conflicting results to a pharmacist or qualified clinician.
28. **ASM DDI Source Hierarchy and Discordance Review**: Evaluate antiseizure-medication DDIs using Lexicomp or another authoritative reference as the primary comparator, secondary references such as Drugs.com as corroboration, and LLM output as unvalidated; report severity and mechanism agreement separately, test and document stability across bounded iterative prompts, disclose uncertainty and source discordance, and require pharmacist review for every discordant result before clinical use.
29. **ASM DDI Validation, Agreement Scoring, and Escalation**: Validate each antiseizure-medication DDI against authoritative references using a source hierarchy that prioritizes Lexicomp or another validated compendium, treats Drugs.com as a secondary comparator, and treats LLM output as unvalidated; score severity and mechanism agreement separately; use predefined, bounded iterative prompts to test answer sensitivity; and mandatorily escalate uncertain, unsupported, or discordant interactions to a pharmacist or qualified clinician.
30. **Reference-Standard DDI Evaluation**: Evaluate DDI outputs against curated databases using an explicit source hierarchy; measure interaction detection, severity agreement, and unsafe omissions; and test whether findings remain stable across predefined, bounded iterative prompting strategies.
31. **ASM DDI Source-Hierarchy and Recommendation-Stability Gate**: Verify every antiseizure-medication interaction against an authoritative interaction database, using Lexicomp as the primary comparator and Drugs.com as a secondary comparator where available; track whether the identified interaction, severity, or management advice changes across each bounded prompt iteration; assess severity separately from management advice; and escalate uncertain or high-risk combinations to a pharmacist before clinical use.
32. **ASM DDI Benchmarking and Provenance Audit**: Benchmark antiseizure-medication DDI findings against Lexicomp as the reference standard, using database-backed sources ahead of model output; run predefined iterative-prompting sensitivity tests; report interaction-severity concordance; track omissions and false-positive interaction claims; and label each finding explicitly as database-backed or model inference.
33. **Reference-Backed DDI Benchmark Procedure**: Compare generated DDI assessments against authoritative references such as Lexicomp; evaluate interaction detection separately from severity classification and management accuracy; run and record predefined iterative-prompt variants to test answer sensitivity; and require source-backed pharmacist review before clinical use.
34. **Source-Grounded ASM DDI Error and Stability Workflow**: For high-risk antiseizure medication combinations, compare each candidate interaction against authoritative drug references, using Lexicomp as the reference comparator and Drugs.com as a secondary comparison source; normalize severity labels and record agreement or disagreement; classify missed reference interactions as potential false negatives and unsupported model-reported interactions as potential false positives; run a predefined, bounded iterative-prompt sequence and document changes in interaction detection, severity, or recommendation; and require pharmacist review before clinical interpretation or use.
35. **Evidence-Backed LLM DDI Validation Workflow**: Validate each LLM-generated interaction against an authoritative comparator such as Lexicomp; score severity and mechanism concordance separately; run predefined, bounded iterative prompts to test answer stability; verify that every interaction claim cites a supporting source; and mandatorily escalate conflicting, uncertain, or unsupported interactions to a pharmacist or qualified clinician.

## Comparative Benchmark: LLM-Based ASM DDI Checks

LLM-generated DDI findings are screening drafts, not standalone clinical evidence. For antiseizure medication benchmarks, compare each generated answer and any Drugs.com output against authoritative databases such as Lexicomp; document prompts and finite iteration limits; classify severity explicitly; record unresolved disagreements without resolving them by model consensus; and route outputs to pharmacist review before clinical use.

## Evaluation and Prompting

1. Establish the source hierarchy before evaluation: use a validated interaction database such as Lexicomp as the authoritative comparator, use Drugs.com as a secondary comparison source, and treat LLM output as unvalidated.
2. Compare severity and mechanism agreement separately across sources, and record omissions as potential false negatives requiring review.
3. Run a predefined, bounded sequence of iterative prompts and record whether each iteration changes the identified interaction, severity, mechanism, or recommendation.
4. Do not interpret prompt-to-prompt agreement or improvement as validation. Never substitute LLM output for Lexicomp or another validated interaction database.

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

*   PubMed PMID 41994367: "Comparative performance of large language models and Drugs.com versus Lexicomp for antiseizure medication drug-drug interactions: A cross-sectional study with iterative prompting analysis." Explor Res Clin Soc Pharm, 2026 Jun. https://pubmed.ncbi.nlm.nih.gov/41994367/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
