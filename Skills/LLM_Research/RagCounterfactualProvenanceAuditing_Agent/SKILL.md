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
name: 'rag-counterfactual-provenance-auditing'
description: 'Audit whether cited RAG claims respond to retrieved evidence using controlled factual substitutions, repeated generation, and unchanged-context controls, while separating citation support from provenance.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# RAG Counterfactual Provenance Auditing

## Overview

Audit the behavioral dependence of cited claims on retrieved passages by replacing a factual value and measuring how generation responds. Citation support and evidence use are different properties: an entailed citation can agree with an answer generated from memory. This workflow produces controlled sensitivity measurements and uncertainty statements, not definitive labels of internal provenance.

## When to Use This Skill

- Evaluate whether a RAG system follows retrieved facts or persists with an original answer when evidence changes.
- Investigate citations that support a claim but may not have influenced its generation.
- Compare models, prompts, or attribution signals on the same claims and controlled interventions.
- Design a research audit with unchanged-context controls and reproducible claim-level records.

## Core Capabilities

1. **Freeze the evaluation setup.** Record model and endpoint versions, system and user prompts, retrieved passage text and order, citation identifiers, decoding settings, seed support, and tool configuration. Hold retrieval fixed for generation-only auditing; disable fresh retrieval or record any unavoidable access as a confound. Use isolated requests with identical history to prevent one condition from contaminating another.

2. **Extract atomic citable claims.** Split baseline answers into independently assessable factual claims. Map each claim to its cited passages, exact supporting spans, subject, relation, and target value. Define citable-claim eligibility before inspecting intervention outcomes; retain exclusions and their reasons. An unsupported, ambiguous, or non-substitutable claim remains part of coverage reporting but need not enter the intervention denominator.

3. **Build valid factual substitutions.** Replace one target value per experimental condition with a controlled false variant of the same semantic type, preserving grammatical fit, units, formatting, and surrounding meaning. Verify the original value and document why the replacement is false relative to the source. Treat proposed variants as unresolved if falsity or coherence cannot be established. Specify substitution distance, such as numeric displacement or entity similarity, without assuming different distance measures are interchangeable.

4. **Check the entire retrieved context.** Inventory repeated mentions, aliases, summaries, tables, and corroborating passages. For the primary condition, replace all occurrences expressing the same target fact consistently, recording every changed span. Preserve unrelated facts, passage order, and citation identifiers. A cited-passage-only intervention is a separate condition: if other passages retain the original fact, persistence cannot establish memory use. Exclude or separately analyze interventions that create contradictions or require broad rewriting.

5. **Run paired controls and interventions.** Generate repeated unchanged-context outputs and repeated substituted-context outputs with matched settings and, when available, matched seeds. Keep the question unchanged; exclude questions that themselves reveal the target answer or mark them as confounded. Randomize condition order and record raw requests, responses, errors, latency, and actual repeat counts. Seeds help reproducibility but do not guarantee deterministic service behavior.

6. **Score semantic response sensitivity.** Align outputs to the original subject and relation rather than relying on string replacement. Label each target response as original value, substituted value, other value, explicit conflict or rejection, abstention, omitted claim, or unscorable. Track citations and textual entailment separately. Use a fixed rubric and human review of ambiguous cases; any automated judge should be blind to the desired conclusion and its errors documented.

7. **Quantify effects and uncertainty.** Report outcome distributions for each condition, baseline instability, and intervention effects with explicit denominators. Preserve unsuccessful or omitted responses as outcomes rather than silently removing them. Cluster uncertainty estimates by question when multiple claims share a question, and avoid treating repeated generations as independent claims. Report counts and descriptive results when the sample is too small for useful inference.

8. **Interpret attribution conservatively.** Following a substituted value demonstrates sensitivity under that intervention, not proof of exclusive contextual provenance in the original answer. Original-value persistence is consistent with memory use but also with conflict resolution, implausibility detection, ignored passages, or failed interventions. NLI measures support, not causal use; context-versus-no-context distribution contrasts are optional exploratory signals, not validated provenance verdicts.

### Execution sequence and time budget

- **Minutes 0–3:** Confirm inputs and model access, freeze configuration, and select a bounded claim sample. Predeclare eligibility, variant rules, repeat counts, outcome rubric, and comparison metrics. Choose a pilot that fits available latency and call limits; the time target is not a statistical power guarantee.
- **Minutes 3–6:** Extract claims and construct substitutions. Review each edit for type consistency, false-value validity, remaining original evidence, and unintended semantic changes. Reject invalid interventions before examining their generated outcomes.
- **Minutes 6–12:** Execute unchanged and substituted conditions with the declared repeats. If feasible, test multiple independently valid variants per claim and summarize results by variant distance. Keep additional variants from overweighting their source claims.
- **Minutes 12–15:** Score outputs, calculate paired effects, and return the audit report. If access, validation, or latency prevents completion, return an explicitly incomplete pilot with completed counts, missing evidence, and a next step; never fill gaps with simulated measurements.

### Measurement definitions

Let A be the original target value and B the substituted value. For each eligible claim and variant, let U denote unchanged-context runs and S substituted-context runs. Compute each probability as its outcome count divided by all completed runs in that condition, with unscorable outputs retained as a separately reported category; report API failures outside that denominator.

- **Variant adoption:** P(B | S), alongside spontaneous B production P(B | U).
- **Adoption effect:** P(B | S) − P(B | U).
- **Original-value persistence:** P(A | S), alongside baseline fidelity P(A | U).
- **Persistence effect:** P(A | S) − P(A | U).
- **Unchanged-context instability:** 1 − P(A | U), decomposed into the other outcome categories. Treat this as variation relative to the selected baseline claim, not a pure measure of randomness.
- **Coverage:** Eligible claims / all extracted citable claims, plus attempted, completed, invalid, and excluded intervention counts and reasons.

Average variants within a claim before aggregating claims, and disclose weighting across questions. Provide paired uncertainty intervals using resampling at the question level when sufficient independent questions exist. Include a sensitivity analysis for scoring ambiguity or unscorable outputs when these could alter the conclusion. Do not convert persistence rates or effect sizes into a percentage of proven parametric provenance.

### Optional attribution comparison

If comparable token probabilities are available, evaluate a predeclared context-versus-no-context distribution contrast on identical candidate claims, documenting token alignment, normalization, and model settings. If probabilities are unavailable, skip this analysis; do not substitute self-reported confidence. No-context answer accuracy may indicate available knowledge but does not establish whether that knowledge caused the contextual answer.

Compare attribution methods on the same eligible claims and disclose method-specific missingness. Any learned thresholds require separate calibration and held-out evaluation. Report precision and false-positive behavior at the observed prevalence before proposing operational actions; do not automatically suppress claims or relabel citations from an exploratory score.

## Inputs / Outputs

**Required inputs**

- Questions, baseline answers, citation mappings, and the exact retrieved passages presented to the generator.
- An authorized generation interface or saved repeated outputs for all conditions, with model and decoding metadata.
- An evaluation budget, claim-selection rule, substitution-validation evidence, and outcome-scoring rubric.

**Optional inputs**

- Multiple model configurations, variant-distance strata, token probabilities, and independent adjudications.

**Outputs**

- A run manifest containing configuration, sampling decisions, intended and actual repeats, condition ordering, and limitations.
- A claim ledger with question ID, claim ID, original claim, citations, source spans, original and replacement values, changed spans, validity checks, and exclusion reasons.
- A response ledger with condition, variant ID, repetition or seed, raw response, aligned outcome, citation-support assessment, scoring uncertainty, and request status.
- An audit report containing coverage, per-condition outcomes, paired effect sizes, uncertainty, baseline variation, representative successes and failures, and unresolved confounds.

Return these as structured sections or machine-readable records within the authorized output destination. A completed audit must make each aggregate traceable to its claim and response records, distinguish support from sensitivity, and state that behavioral evidence does not uniquely identify internal provenance. If only a design is possible, label it as a protocol rather than an executed audit.

The supplied source summary reports study-specific parametric-provenance estimates of 1.2%–15.7%, and 4.6% and 8.5% for two models at its most conservative substitution distance. These are contextual findings, not expected rates, acceptance thresholds, or results of this workflow. The summary also reports chance-level NLI discrimination and limits to actionable detection; do not infer a universal detector accuracy or reconstruct unprovided experimental details.

## References

- [Cited but Unused: Parametric Provenance in Retrieval-Augmented Generation](http://arxiv.org/abs/2611.17161v1), Tanzeel Abbas, Tehseen Zia, and David Windridge; supplied arXiv finding, version 1, published November 12, 2026. This skill is based on the supplied summary; its repeat controls, reporting schema, and time budget are workflow design choices, not a claimed reproduction of the full paper.
