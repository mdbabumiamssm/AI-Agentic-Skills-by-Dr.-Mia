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
name: 'findendro-llm-backtest-contamination-audit'
description: 'Audit temporal knowledge in financial LLM backtests using verified FinDendro protocols, comparing observed factual responses with cutoff assertions and documenting contamination risks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# FinDendro LLM Backtest Contamination Audit

## Overview

Use FinDendro to investigate whether historical financial outcomes are already accessible to an LLM and whether its knowledge claims match observed behavior. The supplied finding describes 2,226 dated U.S. facts spanning 1999–2026 across macroeconomic series, stock prices, and Federal Reserve decisions. This skill produces a reproducible risk assessment; negative recall, historical-date instructions, and vendor cutoff assertions cannot certify an uncontaminated backtest.

## When to Use This Skill

- Before interpreting an LLM's retrospective financial forecasting results.
- When selecting evaluation dates based on a stated training cutoff.
- When a backtest excludes outcomes the model claims to know.
- When testing whether an instruction to answer as of a historical date prevents hindsight leakage.
- When comparing temporal knowledge across model snapshots or reviewing contamination claims in an experiment.

## Core Capabilities

1. **Verify the complete protocol before execution.** Read the versioned paper, methods, appendices, and any officially linked benchmark artifacts. Record source version, retrieval date, artifact identifiers, and relevant section references. Verify fact construction, date semantics, sampling, candidate generation, option order, all prompt templates, controls, scoring denominators, uncertainty estimation, and boundary rules. Do not reconstruct a purported exact replication from the abstract alone. If necessary artifacts or details are unavailable, return a protocol-gap assessment and an explicitly provisional plan; do not invent data or report benchmark completion.

2. **Scope the experiment and freeze its configuration.** Record the model provider, exact snapshot, query date, decoding settings, system prompt, claimed cutoff with supporting source, and proposed backtest interval. Distinguish vendor claims, model self-reports, release dates, and behavioral estimates. Define the tested domains and date coverage before seeing results. Treat a moving model alias as a reproducibility limitation.

3. **Check dated facts and information availability.** Validate benchmark records against their documented sources. Separate reference periods from first public availability and later revisions. Inspect units, rounding, missing values, and price-adjustment conventions. Keep ground truth and scoring metadata outside the tested model's context. Independently inspect the backtest's feature timestamps and retrieval corpus; a factual-knowledge audit does not replace those checks.

4. **Run aligned behavioral probes.** Use the verified four-prompt design: value-only, undated Recall with knowledge claims, undated Predict with knowledge claims, and historically dated Predict. Preserve shared candidate values and the paper's option construction. Run prompts independently, without browsing, retrieval, or conversational carryover in the tested model. Save exact prompts, option mappings, raw responses, parsing outcomes, and failures. Auditor access to sources is separate from model access during probing.

5. **Apply controls and compare claims with answers.** Follow the published persistence exclusion for numerical facts and class balancing for policy outcomes. Keep controlled and full-sample denominators explicit. Measure factual accuracy, claim frequency, conditional accuracy for claimed knowledge and claimed ignorance, and paired changes under historical dating using the verified definitions. Do not substitute claim-bearing answers for a separately defined value-task score. Recompute the appropriate chance reference for the actual sample; the finding's approximately 0.27 baseline is not universal. Report invalid responses and exclusions rather than silently dropping them.

6. **Estimate temporal boundaries with uncertainty.** Implement the paper's verified aggregation, smoothing, thresholds, baseline period, and resampling procedure, with a source pointer for every methodological choice. Show sample counts and uncertainty by time period and domain. If the historical baseline is near chance, coverage is sparse, or observations end before a decline, mark the boundary unidentifiable or censored as appropriate. A single incorrect response or nonsignificant period does not establish where knowledge ends. Label any alternative estimator or subset analysis as an adaptation.

7. **Translate evidence into a backtest risk assessment.** Compare the proposed interval with measured knowledge, knowledge claims, and asserted cutoffs. Flag overlap with demonstrated factual knowledge as contamination risk, without claiming to identify particular training examples or quantify inflated trading returns. Describe a negative test as absence of detected knowledge under this protocol. Recommend prospective evaluation, independently documented point-in-time model snapshots, or further targeted probes where warranted; any later candidate window remains conditional on coverage and uncertainty.

8. **Deliver a bounded, reviewable result.** Within the time budget, analyze existing responses or execute a feasible, prespecified subset. Clearly distinguish full replication, partial audit, and protocol-gap assessment. Do not label a pilot as the full benchmark or imply that a time limit validates incomplete evidence. Identify unfinished work and the exact missing inputs needed to continue.

## Inputs / Outputs

**Inputs**

- Versioned paper and appendices; available official benchmark data, question files, and implementation.
- Model snapshot or existing raw response logs, configuration, and vendor cutoff evidence.
- Proposed backtest dates, target outcomes, forecast origins, feature availability, and retrieval/tool settings.
- Execution budget and any approved sampling constraints; a full live run is not required to review existing evidence.

**Outputs**

- **Protocol manifest:** verified sources and versions, methodological settings, unresolved details, and deviations from the paper.
- **Response ledger:** fact identifier, domain, reference date, public availability date where available, model snapshot, prompt variant, declared date, option mapping, raw response, parsed value and claim, correctness, control eligibility, and error status.
- **Results table:** period/domain counts, controlled factual scores, claim rates, conditional scores, dated-versus-undated comparisons, baseline definitions, uncertainty, and missingness. State each numerator and denominator and distinguish paper-reported results from newly measured results.
- **Boundary comparison:** vendor assertion, model assertion, behavioral estimate, claim-based estimate, coverage limits, and uncertainty; use “not estimable” where justified.
- **Risk memo:** audited scope, evidence of knowledge overlapping the backtest, unresolved contamination channels, proposed evaluation periods with rationale, and next steps. Include: “Failure to detect recall does not prove that this backtest is uncontaminated.”

An output is valid when its scope and provenance are explicit, reported measurements are traceable to actual responses, and limitations prevent unsupported certification. If execution is blocked, provide the protocol gaps and unexecuted plan without fabricated scores. Keep outputs in the response unless the user requests saved artifacts.

## References

- [LLMs Don't Know What They Know: Benchmarking Temporal Self-Knowledge for Financial Backtesting, arXiv:2611.06808v1](http://arxiv.org/abs/2611.06808v1) — source finding; use this pinned version for protocol verification.
- [Versioned full text and appendices](https://arxiv.org/html/2611.06808v1) — inspect Sections 3–5 and the appendices before implementing the audit. Follow repository or dataset links only when their association with this paper is verified.
