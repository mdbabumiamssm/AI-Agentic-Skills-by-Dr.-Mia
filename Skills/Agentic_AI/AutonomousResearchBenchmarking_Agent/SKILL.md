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
name: 'autonomous-research-benchmarking'
description: 'Set up and assess autonomous language-model research using Auto-Research Arena principles, separating admissibility, target improvement, resource budgets, and independently validated scientific progress.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Autonomous Research Benchmarking

## Overview

This skill prepares reproducible evaluations of autonomous research systems in language modeling and assesses their experimental evidence. Auto-Research Arena separates quality and resource requirements for admissibility from the target used to measure improvement, with nine targets spanning model, data, training, and inference efficiency ([paper](http://arxiv.org/abs/2610.35866v1)). The workflow below applies that distinction while keeping benchmark success separate from independently validated scientific progress.

## When to Use This Skill

- Configure an Auto-Research Arena evaluation or audit an existing experiment bundle.
- Compare research harnesses under a shared environment and explicit search budgets.
- Evaluate efficiency improvements while enforcing model-quality requirements.
- Assess whether a discovered idea warrants replication, ablation, or transfer testing.
- Produce a reproducible experiment plan when compute or release artifacts are unavailable.

## Core Capabilities

1. **Verify the benchmark artifacts.** Read the versioned paper and the linked official repository. Inspect available documentation, target configurations, evaluation code, and release artifacts before deriving commands or thresholds. Pin a release and resolved commit when available, otherwise a specific commit; record retrieval date and relevant artifact hashes. A reachable repository alone does not establish that a complete executable release exists. Mark missing artifacts as unverified and produce a bounded plan rather than inventing interfaces or claiming an official run.

2. **Define the evaluation contract.** Select the requested target from the verified configuration and record its metric, units, optimization direction, baseline, quality thresholds, and resource limits. Specify allowed modifications and protected evaluation components. Keep each target separate unless the release defines an aggregation method; do not average incompatible units or invent the nine target definitions from the abstract. Label user-defined extensions explicitly.

3. **Separate qualification from improvement.** Evaluate every candidate against all required quality and resource constraints before ranking its target score. Use `admissible`, `inadmissible`, or `unverified`, with evidence and reasons. Missing measurements are unverified, and an attractive score cannot compensate for a failed constraint. Report improvement against an admissible baseline using the official scoring rule; if none exists, report raw values and signed differences with a stated direction, without presenting these as official scores.

4. **Control the comparison and budgets.** Freeze the data splits, preprocessing, tokenizer, evaluator, environment, and target configuration. Record the research foundation model and settings separately from the language model being trained. Compare harnesses with matched access and predeclared budgets for attempts, wall time, accelerator time, model calls or tokens, and cost where applicable. Distinguish total search expenditure from candidate training and inference resources. Include failed attempts and retries in accounting, document any official exclusions, and stop at the first binding limit. Include a simple sequential search baseline where feasible; do not infer superiority from harness complexity.

5. **Run or audit reproducibly.** Use commands verified against the pinned documentation, first checking a baseline or small smoke run. Record seeds, hardware, dependency versions, data identifiers or hashes, prompts, configuration, code changes, timestamps, logs, exit status, resource consumption, and checkpoint or artifact identifiers. Preserve the proposal-to-experiment history, including rejected ideas and failures. Keep held-out evaluation data out of agent feedback according to the benchmark protocol, and evaluate submitted artifacts with the protected evaluator. Label modified or reduced runs as smoke tests or custom evaluations.

6. **Assess outcomes and scientific evidence.** Report admissibility rate, best admissible target result, baseline difference, total search resources, and failure categories for each target and harness. For repeated runs, give seed counts and variation; identify single-run findings as preliminary. Freeze selected discoveries before independent reruns, use fresh seeds and held-out evaluation, and test component ablations and transfer to another suitable setting. For quality claims, state whether data and FLOPs are matched for training-loss comparisons and whether model size is matched for held-out bits-per-byte (BPB) comparisons. Report remaining differences. Benchmark improvement supports a claim within its protocol; broader scientific progress requires independent supporting evidence.

7. **Finish a bounded assessment within 15 minutes.** Use this time for artifact inspection, contract preparation, analysis of supplied results, and a smoke test only when feasible. Full research searches and replication may require separate execution budgets. Return a valid assessment marked `plan-only`, `smoke-tested`, `evaluated`, or `blocked`, with actual evidence, missing prerequisites, and the next executable step. Never label an unfinished training run or unexecuted plan as a completed benchmark evaluation.

## Inputs / Outputs

**Inputs**

- Requested target or research objective and whether the task is setup, execution, or audit.
- Official artifact location and version, or permission to retrieve publicly available artifacts.
- Harness definitions, foundation-model settings, baseline, dataset specifications, and evaluation protocol.
- Hardware and environment details, authorized experiment budgets, and stopping conditions.
- Existing configurations, logs, checkpoints, metrics, and experiment history when available.

If essential inputs are missing, identify them explicitly and complete the portions supported by available evidence. Do not silently supply benchmark thresholds, claim unavailable compute, or start a full search without a defined budget.

**Outputs**

- **Evaluation contract:** source versions, target definition, admissibility rules, allowed changes, baseline, comparison controls, budgets, and stop conditions.
- **Reproduction manifest:** exact commands when verified, configuration and environment identifiers, data provenance, seeds, and pointers to evidence artifacts. Omit credentials from logs and reports.
- **Results ledger:** one row per attempt with target, harness, seed, status, quality measurement, resource measurements, target score, admissibility decision, failure reason, and evidence pointer. Use explicit missing values rather than fabricated measurements.
- **Assessment:** per-target findings, search costs, uncertainty, comparability limits, execution status, and a clear distinction between observed results and proposed work.
- **Validation plan:** independent reruns, ablations, held-out checks, and matched-resource comparisons needed to substantiate any proposed scientific claim.

The output is valid when every reported result has traceable evidence, every ranked result satisfies the declared constraints, and incomplete or incomparable results are visibly labeled. For plan-only work, provide the complete proposed contract and outstanding dependencies without measured-performance claims.

## References

- [Auto-Research Arena: Benchmarking and Advancing Auto-Research in Language Modeling, arXiv:2610.35866v1](http://arxiv.org/abs/2610.35866v1) — source finding and benchmark framing; submitted October 26, 2026.
- [Auto-Research Arena official repository](https://github.com/Auto-Research-Arena/auto-research-arena) — paper-linked source for verifying available code, configurations, evaluation rules, and release artifacts; pin the version used in each evaluation.
