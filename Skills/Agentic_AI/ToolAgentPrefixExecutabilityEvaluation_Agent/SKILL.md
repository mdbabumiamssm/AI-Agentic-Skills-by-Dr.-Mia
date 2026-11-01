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
name: 'tool-agent-prefix-executability-evaluation'
description: 'Evaluate tool-agent trajectories against schemas, state preconditions, and data dependencies, reporting executable-prefix reliability and recovery alongside final success.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Tool-Agent Prefix Executability Evaluation

## Overview

Evaluate intermediate tool-call feasibility alongside terminal task success. Check schemas, state preconditions, and data dependencies at each prefix, then report PEC, FED, TER, normalized AUC-PEC, and recovery to expose reliability differences hidden by outcome-only evaluation.

## When to Use This Skill

- Audit successful episodes containing rejected calls or repaired arguments.
- Compare tool agents, planning strategies, or pre-execution validators.
- Diagnose the earliest point where a workflow becomes non-executable.
- Complement computer-use completion evidence with tool-trajectory checks.

## Core Capabilities

1. **Evidence preparation.** Inventory trace completeness, schema versions, initial snapshots, state transitions, and outcome labels. Record task IDs and evaluator versions in a manifest; preserve the original logs.
2. **Deterministic validation.** Bind each call to its actual tool interface; check required fields, types, validators, reference provenance, object state, and prerequisites. Use an isolated replay environment or authoritative recorded transitions. Never replay mutations against production merely to evaluate them.
3. **Prefix accounting.** Identify the first rejection and retain evidence for its failed predicate. Preserve the complete observed episode for outcome and recovery analysis.
4. **Metric computation.** Apply the definitions below with explicit counts, horizons, and denominators.
5. **Comparative diagnosis.** Pair task IDs across methods, inspect failure examples, and distinguish reliability changes from changes in trajectory length or final success.

### Workflow

- **Minutes 0–3:** Select a tractable existing trace cohort and record coverage. Resolve missing schemas or snapshots before scoring; label unverifiable episodes unresolved rather than silently passing or failing them.
- **Minutes 3–9:** Build an action ledger containing tool, arguments, available bindings, pre-state, validation evidence, accepted output, and post-state. Check references against information available at that step, never future observations. Advance replay state only after accepted transitions. Stop first-failure replay at rejection; retain later observed events separately.
- **Minutes 9–12:** Compute metrics on the declared evaluable cohort. Report unresolved and excluded counts outside its denominator. If final outcomes are missing, report RR as unavailable rather than substituting failure labels.
- **Minutes 12–15:** Check arithmetic invariants and produce the report. If the workload exceeds the time budget, return an explicitly partial cohort report and remaining work; do not claim complete evaluation.

### Metric contract

Use Section 2.2 of the reference. For N rollouts, T_i counts emitted tool-action slots; Y_i is benchmark success. Exclude planner/verifier messages. Represent no-tool termination by one slot, valid only when benchmark semantics allow it.

Let F_i = FED_i be the first invalid depth, or T_i + 1 when all calls pass. With I denoting an indicator:

| Output | Empirical definition |
| --- | --- |
| PEC(k) | sum_i I(F_i > k) / N |
| FED | Each F_i and mean_i F_i |
| TER | sum_i I(F_i = T_i + 1) / N |
| Normalized AUC-PEC | mean_i ((F_i - 1) / T_i) |
| RR-Eligible | B = sum_i I(F_i <= T_i) |
| RR | sum_i Y_i I(F_i <= T_i) / B; undefined if B = 0 |

Report n_k = sum_i I(T_i >= k). Optional eligible-depth PEC divides the PEC numerator by n_k; label it separately. Short completed trajectories contribute zero beyond their horizons. Later recovery never resets FED. AUC-PEC is not trapezoidal integration of the aggregate curve.

Use first-failure precedence: schema, dependency, precondition, state-consistency, information (Appendix C). Keep secondary evidence separately.

### Verification and interpretation

- Check that PEC never increases, 0 <= TER <= AUC-PEC <= 1, and B = N(1 - TER), allowing rounding tolerance.
- Include a small hand-checked fixture covering first-call failure, later failure, error-free completion, and successful recovery. Use synthetic fixtures only for implementation checks, never as benchmark results.
- Cross-tabulate final success against full executability. Report the recovery numerator and denominator, horizon distribution, and first-failure evidence so readers can audit the result.
- For comparisons, hold tool versions and execution budgets constant. Report paired task-level uncertainty when feasible; do not infer significance from overlapping or separated marginal intervals alone.
- Keep plan validation decisions separate from observed execution evidence. Record validation cost when comparing an intervention, and avoid claiming that validator acceptance proves executability.

## Inputs / Outputs

**Inputs:** Ordered episode logs with task/variant IDs, tool calls and responses, tool schemas and validators, initial state and state evidence, reference bindings, termination reasons, benchmark outcome labels, and execution limits. Optional inputs include plan-validation logs and matched comparison cohorts.

**Outputs:** A cohort manifest; per-action evidence ledger; per-episode T_i, FED_i, outcome, and first-failure cause; a depth table containing PEC and n_k; aggregate TER, mean FED, normalized AUC-PEC, RR with counts, and final success; plus coverage limitations and diagnostic examples. Return tables inline unless the user requests saved artifacts.

## References

- [Prefix Executability: Evaluating Tool-Using Agents Beyond Final Success, arXiv:2610.40444v1](http://arxiv.org/abs/2610.40444v1). Metric definitions: Section 2.2; replay and classification: Appendix C; validation accounting: Appendix D. [Full text](https://arxiv.org/html/2610.40444v1).
