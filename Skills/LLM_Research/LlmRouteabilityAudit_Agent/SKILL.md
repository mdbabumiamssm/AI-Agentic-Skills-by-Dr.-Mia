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
name: 'llm-routeability-audit'
description: 'Audit existing benchmark-by-model scores for routing headroom, tied best outcomes, and gains over the single best model before interpreting LLM router comparisons.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# LLM Routeability Audit

## Overview

Audit an evaluated LLM portfolio for prompt-specific opportunities before investing in routing or interpreting router comparisons. Following the retrospective scope of RAP and RouteabilityBench, distinguish available score headroom from demonstrated, learnable routing gains while exposing ties and portfolio dependence. This audit requires an existing score matrix and cannot certify an unevaluated portfolio.

## When to Use This Skill

- A benchmark-by-model score matrix exists and the team must decide whether router development is justified.
- Coarse scoring produces many tied best models or the strongest fixed model appears difficult to improve upon.
- A routing experiment claims gains over a single-model baseline.
- Portfolio composition, benchmark selection, or predictor choices may explain apparently inconsistent routing results.
- A report needs routeability diagnostics and statistical evidence alongside routing accuracy.

## Core Capabilities

Execute the following workflow within the 15-minute audit budget using existing results. The computational definitions and reporting defaults below operationalize the supplied finding; they are not a claim to reproduce every detail of the published RAP implementation.

1. **Validate the evaluation matrix.** Read prompt IDs, model IDs, scores, scoring rules, and portfolio metadata. Require comparable scores on the same prompts; record metric direction, score granularity, evaluation settings, missingness, and duplicate IDs. Never silently treat missing evaluations as incorrect answers. Use a documented common evaluated subset or return an incomplete audit if comparable data cannot be established. Analyze benchmarks separately unless an aggregation rule was specified in advance.

2. **Compute fixed-model performance and oracle headroom.** For a higher-is-better matrix `S[i,j]` with `n` prompts, calculate each model's mean score, the descriptive single best model score `B = max_j mean_i S[i,j]`, and the per-prompt oracle score `O = mean_i max_j S[i,j]`. Report absolute headroom `H = O - B` in score units. For lower-is-better metrics, explicitly transform the objective first. The oracle uses observed outcomes and is an upper bound on this matrix, not a deployable router or proof of predictability.

3. **Audit tied best outcomes.** Define `T_i` as all models achieving the best score for prompt `i`; use exact equality for discrete scores and a documented tolerance for continuous scores. Report the distribution of `|T_i|`, fractions of unique-best and multiple-best prompts, and all-model ties. Split all-model ties by score level so universal failure is not confused with universal success. Report each model's unique-best frequency and tied-best membership frequency; the latter need not sum to one. Do not create arbitrary winner labels to manufacture apparent routing diversity.

4. **Locate actionable opportunities.** For a recorded choice among equally strong fixed models, compute each prompt's oracle improvement over that model. Report the fraction with positive improvement and its magnitude, along with how frequently the fixed model belongs to `T_i`. Inspect whether gains are concentrated in a small subset or depend on portfolio membership. Document inclusion criteria, model versions, family composition, and any subsets tried; outcome-driven portfolio selection must be labeled exploratory and validated independently.

5. **Assess routers without evaluation leakage.** If router selections are supplied, compute their realized scores on the same held-out prompts and paired differences from the baseline. Select a deployment baseline on training or validation data and freeze it before test evaluation; report the test-matrix descriptive best model separately. Keep test outcomes out of router training, feature tuning, and portfolio selection. Report tied-best hit rate as a diagnostic alongside achieved score, never as a substitute for score gains. Without held-out router outputs, limit conclusions to retrospective opportunity.

6. **Quantify uncertainty and multiple comparisons.** Report effect sizes and paired uncertainty intervals, resampling at the independent unit when prompts share groups. For confirmatory router comparisons, use a test suited to the score and dependence structure, such as a paired permutation test where exchangeability holds. Declare the hypothesis family across the comparisons supporting the claim and apply a family-wise correction such as Holm; retain raw and adjusted p-values, significance level, and the full comparison count. Distinguish multiplicity correction from model-family dependence. If splits, replication, or independence assumptions are missing, state that inferential claims remain unsupported.

7. **Issue an evidence-bounded gate.** Use a task-specific minimum useful gain supplied before inspecting results. If none is provided, report headroom without inventing a universal threshold. Return one of: `insufficient data`, `no useful observed headroom`, `opportunity present; predictability unproven`, or `held-out routing gain supported`. Reserve the final status for a practically relevant held-out gain with appropriate corrected evidence. Scope every status to the evaluated benchmark, scoring procedure, portfolio, and split; do not infer that adding embeddings, clustering, regressors, or abstention will reliably improve routing.

## Inputs / Outputs

**Required inputs**

- A score matrix or long-form table containing `prompt_id`, `model_id`, and numeric `score`, with benchmark identity where applicable.
- Metric definition, score direction, granularity, missing-value policy, and evaluation configuration.
- Portfolio membership and construction rationale, including model versions when available.

**Optional inputs**

- Train/validation/test assignments, repeated-run identifiers, and prompt grouping information.
- Held-out router choices keyed by prompt, baseline selection procedure, and all compared router configurations.
- A predefined minimum useful gain, significance level, and comparison family.
- An existing RouteabilityBench checkout or verified package location for exact protocol reproduction.

**Outputs**

- A concise audit report containing data coverage, assumptions, fixed-model scores, oracle score, headroom, tie diagnostics, portfolio construction, and the gate status with rationale.
- A machine-readable diagnostic table containing per-prompt best score, tied-best model IDs and count, baseline score, and available improvement; add router scores and paired gains when provided.
- A statistical summary with effect sizes, uncertainty intervals, correction method, hypothesis family, and unresolved limitations, or an explicit statement that router inference was not performed.
- A reproducibility record containing input identifiers or hashes, split definitions, tie tolerance, aggregation rules, seeds, and commands actually executed.

Use `read_file` for local results and `run_shell_command` for bounded calculations. Use `web_fetch` to verify the source and any official implementation link before invoking package-specific commands; do not invent a repository URL or RAP command-line interface. Exact replication requires verified protocol details and an implementation version; otherwise label the result as a finding-grounded audit. If inputs or computation exceed the time budget, return an incomplete report with missing checks rather than fabricate a completed gate.

Before delivering, check that `O >= B` within numerical tolerance, every prompt has at least one tied-best model, and all reported fractions use explicit denominators. Preserve original score units and distinguish descriptive diagnostics from confirmatory results.

## References

- [Routeability Before Routing: Routeability Audit Protocol (RAP) and RouteabilityBench for Audited LLM Model Selection](http://arxiv.org/abs/2611.34825v1), arXiv:2611.34825v1, published November 25, 2026. Source for the retrospective audit scope, tied-best analysis, single-best-model comparison, portfolio dependence, and family-corrected reporting requirements.
