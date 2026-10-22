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
name: 'quantitative-research-agent-evaluation'
description: 'Evaluate quantitative research agents with Quant-Arena through executable signals, research-only feedback, frozen strategies, and cost-aware held-out A-share assessment.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Quantitative Research Agent Evaluation

## Overview

Use Quant-Arena to assess an agent's complete cycle of market hypothesis, executable signal construction, historical testing, and strategy selection. The supplied finding describes research feedback restricted to 2022–2023 and one frozen strategy assessed on sealed 2024 A-share data, with trading costs and constraints applied by a cash-account evaluator. This skill supports reproducible evaluation and separates leaderboard performance from evidence of out-of-period profitability.

## When to Use This Skill

- Set up a Quant-Arena evaluation for a model or agent harness.
- Prepare feature-generation and stock-scoring code for the benchmark's submission interface.
- Check research versus held-out separation, data leakage, and execution assumptions.
- Compare agents under consistent research budgets, data access, and evaluator versions.
- Interpret held-out results, including whether a strategy exceeds the benchmark market return.

## Core Capabilities

1. **Establish the benchmark contract.** Retrieve the versioned paper in References and follow its official project links to available code, data instructions, and submission documentation. Record the evaluator revision, data release, interface, dependency requirements, scoring definitions, and permitted agent tools. Do not invent installation commands, repository addresses, score weights, or API signatures. If executable resources are unavailable, provide a readiness report identifying the missing prerequisites instead of claiming a completed benchmark run.

2. **Prepare a reproducible evaluation.** Use the documented environment and official example strategy to check loading, feature computation, stock scoring, and evaluator invocation. Record the model identifier, harness version, prompts, random seeds where supported, compute limits, and research-call budget. Keep these conditions consistent across comparisons and document any deviations from the published protocol.

3. **Enforce temporal separation.** Permit strategy-development feedback only from 2022–2023; reserve 2024 evaluation for the frozen submission. Check that features use information available at the scoring timestamp, transformations are fitted only on permitted data, and labels or future returns never enter signal inputs. Check timestamp joins, corporate-action adjustments, universe construction, missing-data handling, and survivorship assumptions against the actual data contract. These are audit checks, not assertions that the released benchmark implements particular controls. Isolate sealed data and its feedback from the research agent; record any prior exposure that compromises an independent holdout claim.

4. **Run the research loop.** State a testable market hypothesis and its expected signal direction, then implement features and stock-scoring code using the verified interface. Validate output shape, identifiers, timestamps, finite values, and deterministic behavior where required by that interface. Test on the research period, log the candidate and permitted feedback, and decide whether to revise or retain it. Track all trials and selection decisions so repeated testing is visible; do not select candidates using 2024 outcomes.

5. **Audit trading realism.** Inspect the evaluator's cash accounting, order timing, position sizing, tradability rules, and cost configuration before interpreting results. Verify fees, slippage, and the simulated-order-book market-impact mechanism against the released implementation. Check whether and how A-share constraints such as lot sizes, price limits, settlement restrictions, and suspensions are represented; do not silently substitute assumed rules. Report unresolved details and distinguish official settings from optional sensitivity analyses. Evaluate the reference market using the benchmark's documented convention.

6. **Freeze and evaluate one strategy.** Select the final candidate using research-period evidence, then record a content hash or immutable revision for code, parameters, feature definitions, dependencies, and configuration. Submit it through the documented evaluator interface without permitting code changes or training from held-out feedback. Keep private evaluation credentials and sealed data outside the research agent's access. If results trigger a revision, label that work exploratory and require a new untouched period for an independent assessment; do not present retuning on 2024 as sealed evaluation.

7. **Report performance with the correct scope.** Report the official aggregate score and its available components using verified definitions. The finding describes realized excess return as dominant, with signal quality and feature diversity also contributing; do not infer their numerical weights. Separate research and held-out results, show strategy and reference returns, and identify costs and failed runs. Add risk or turnover diagnostics only when computable and label them as supplementary. A higher aggregate score does not establish market outperformance: the supplied finding reports that none of the evaluated models beat the equal-weight market in 2024. Limit conclusions to the tested A-share universe, period, and execution assumptions.

## Inputs / Outputs

**Inputs**

- Versioned Quant-Arena source and accessible official evaluator, data, and interface documentation.
- Research-period data for 2022–2023 and a controlled mechanism for sealed 2024 assessment.
- Model or agent harness, permitted tools, resource budget, and reproducibility settings.
- Hypotheses, candidate feature and stock-scoring implementations, and any existing experiment logs.
- Evaluator configuration, including reference-market definition, constraints, and transaction costs.

**Outputs**

- A setup manifest with verified versions, commands, access boundaries, and unresolved requirements.
- A research experiment log containing hypotheses, candidate revisions, research-only feedback, and selection rationale.
- One frozen executable strategy with parameters, feature definitions, and an immutable identifier when prerequisites permit execution.
- An evaluation report separating official scores, realized returns, reference returns, supplementary diagnostics, and leakage or execution limitations.
- A clear status: completed evaluation, partial validation, or blocked setup, with evidence and the next required action.

**Time-bounded workflow and acceptance**

Target 15 minutes for a bounded invocation: first verify resources and protocol, then validate an executable candidate and run the available authorized evaluation, reserving time to report results. Full research campaigns or sealed scoring may exceed this window; record partial progress or a blocker explicitly rather than fabricating results or weakening the holdout boundary. A completed evaluation requires successful interface validation, documented research-only selection, an immutable final strategy, and evaluator-produced results traceable to that strategy. A readiness report is a valid bounded output but is not a completed benchmark evaluation.

## References

- [Quant-Arena: Benchmarking LLM Agents for Quantitative Trading Research, arXiv:2610.30113v1](http://arxiv.org/abs/2610.30113v1), submitted October 21, 2026. Source for the benchmark scope, research/held-out split, evaluator overview, scoring components, and reported out-of-period finding; consult the paper and its official project links for executable details.
