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
name: 'strategist-inference-strategy-evolution'
description: 'Evolve reusable executable inference blueprints through development-set evaluation, library reuse, compute budgets, and held-out validation when improving agent reasoning workflows.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Strategist Inference Strategy Evolution

## Overview

Improve an agent's inference procedure through reusable code, measured experiments, and accumulated experience while keeping its models frozen. Strategist retrieves blueprints, proposes candidates, evaluates them, and accumulates reusable structures across tasks; this skill turns that approach into a bounded development workflow with explicit evidence and leakage controls.

## When to Use This Skill

- Replace a manually designed reasoning procedure with an evaluated alternative.
- Reuse successful inference structures for a related task or another model backend.
- Compare routing, verification, revision, or aggregation strategies under a fixed compute allowance.
- Build a persistent strategy library that records failures as well as successes.
- Diagnose whether strategy-search gains survive held-out evaluation.

## Core Capabilities

1. **Establish artifact fidelity.** Read the versioned paper, especially Sections 2–3 and Appendix A, before implementing its interfaces. The paper describes Python blueprints exposing `execute(llm, task) -> str`, with task-specific prompts supplied separately. Its library enforces a common interface, backward-only dependencies, acyclicity, and deduplication. The linked code repository returned HTTP 404 during authoring; do not claim to have verified its implementation. If code becomes accessible, inspect and pin a revision before using repository-specific commands or APIs. Otherwise implement only paper-supported behavior and label adapters as local conventions.

2. **Define a bounded experiment.** Before any model call, record the task, output schema, scorer, baseline, model configuration, sample identifiers, seed, and selection rule. Set ceilings for elapsed time, candidates, model calls, tokens, monetary cost, and per-sample runtime. Reserve evaluation time and cost before spending on search. A practical local allocation is 2 minutes for setup, 8 for exploration, 3 for final evaluation, and 2 for reporting; these timings are not paper settings. Stop launching work when the remaining allowance cannot cover it. Missing credentials, data, or evaluator support must produce a clearly marked blocked or incomplete report, never synthetic measurements.

3. **Separate feedback from final evidence.** The paper distinguishes training traces, aggregate validation feedback, a hidden selection split, and an untouched test split. Implement these as separate access paths: the designer sees development examples and diagnostics, the selector receives hidden aggregate measurements only after proposals are frozen, and the final evaluator owns test labels. Split by task-relevant groups and remove duplicates across partitions. Hidden selection data are still selection data, not unbiased test evidence. Exclude test content and labels from prompts, retrieval indexes, caches, and experience notes. When data are insufficient, report exploratory results and the missing validation rather than claiming generalization.

4. **Retrieve and compose blueprints.** Search existing structures by task requirements, failure modes, and measured cost as well as similarity. Include a structurally different candidate when budget permits. The paper distinguishes new, composed, and reused blueprints; its library retains structures and measurement notes rather than task prompts. For a local registry, assign immutable IDs, source hashes, dependency IDs, interface versions, and provenance. Reject cycles and forward references, deduplicate equivalent structures, and retain each evaluation record even when it refers to an existing blueprint. Keep task prompts in the experiment record, outside the reusable blueprint registry.

5. **Execute through controlled adapters.** Implement the verified entrypoint and return a normalized answer string; collect metrics in the runner instead of changing the interface. Define the local `llm` adapter and task schema explicitly. Keep evaluator labels out of the strategy's task input. Map logical routing tiers to configured backends and record unavailable capabilities rather than silently substituting models. Inspect generated code before execution and run it with timeouts, limited filesystem access, and only the tools and network access required by the task. Propagate one shared budget through nested blueprint calls, retries, and tool calls so composition cannot bypass limits.

6. **Evaluate, revise, and account for compute.** Run the baseline and candidates on identical development samples with the same scorer, answer extraction, and tool permissions. Check schema validity and termination before larger evaluations. Record sample-level correctness or task score, execution failures, calls, tokens, latency, and cost; count failures in the denominator. Separate designer/search expense from deployed inference expense, including retries and failed runs. Use actual configured prices or a declared compute proxy, and never label a proxy as dollars. Revise from development diagnostics and validation aggregates, then archive every attempt with its parent, change rationale, results, and limitations. Do not present unmeasured candidates as successful library entries.

7. **Select and validate without moving the goalposts.** A suitable local rule is to maximize the declared quality metric subject to the cost ceiling, with a predeclared quality tolerance and cost tie-breaker. Freeze the candidate pool and rule before hidden selection. Keep the baseline if no candidate qualifies. Run the selected strategy and baseline once on the untouched test split, reporting sample counts, paired differences, failures, cost, and uncertainty where feasible. Do not revise or reselect using those results; further optimization requires a fresh final holdout. If time or budget prevents this step, mark the candidate as unvalidated. For transfer tests, freeze the structure and declare any adapter or prompt changes; task-specific search disqualifies a claim of transfer without further search.

## Inputs / Outputs

**Inputs**

- Task specification, permitted tools, required answer format, and executable scorer.
- Development data plus distinct selection and final test partitions, or permission to derive disjoint partitions from sufficient data.
- Baseline implementation, configured model access, pricing or compute proxy, and all execution budgets.
- Optional blueprint library with source code, dependency metadata, measurement history, and provenance.

**Outputs**

- Selected executable blueprint or retained baseline, its dependency closure, and separately recorded task prompts and adapter configuration.
- Candidate ledger covering successful, failed, deduplicated, and unevaluated attempts; reusable library updates within the user's requested scope.
- Reproducibility record containing source revision or hash, split identifiers, scorer version, model settings, selection rule, and stopping reason.
- Evaluation report separating development, selection, and final test results; inference and search costs; limitations; and a status of validated, exploratory, incomplete, or blocked.

**Completion check:** Within the 15-minute workflow budget, return the artifacts that were actually produced and their status. Claim a validated improvement only when executable output, comparable baseline measurements, and untouched-test evidence support it. The time target is a reporting deadline, not a guarantee of better accuracy or lower cost.

## References

- [Strategist: Evolving Reasoning with Agents, arXiv:2611.03618v1](http://arxiv.org/abs/2611.03618v1) — source finding.
- [Versioned full text](https://arxiv.org/html/2611.03618v1) — blueprint interface, search loop, split protocol, and implementation appendices. Operational budget defaults and additional safeguards above are local workflow recommendations.
- [Paper-linked Strategist repository](https://github.com/zlab-princeton/Strategist) — linked by the full text; returned HTTP 404 when checked on 2026-11-03. Availability and implementation remain unverified.
