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
name: 'workflowarena-orchestration-evaluation'
description: 'Prepare and evaluate WorkflowArena runs to compare executable multi-agent orchestration through verified deliverables, delegation traces, and reproducible model-harness configurations.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# WorkflowArena Orchestration Evaluation

## Overview

Evaluate how agents construct executable workflows and integrate delegated work using WorkflowArena, a benchmark of 145 long-horizon tasks across eight professional domains. Connect hidden-verifier outcomes with execution traces to distinguish deliverable quality from delegation behavior. The 15-minute outcome covers readiness assessment or analysis of existing results; a full benchmark execution is a separate, longer operation.

## When to Use This Skill

- Prepare WorkflowArena task packages and identify missing execution dependencies.
- Compare model, harness, or effort configurations on identical tasks.
- Diagnose coordination failures in otherwise plausible component outputs.
- Audit delegation counts, verifier results, and resource accounting.
- Design controlled comparisons of workflow sources with a fixed executor.

## Core Capabilities

1. **Acquire and verify the released benchmark.** Read the [official repository](https://github.com/Einsia/WorkflowArena) at a pinned commit. The release provides data, evaluators, and aggregate results; execution frameworks and individual experimental records are excluded. Record this limitation before promising reproduction. Use Python 3.11+, Git, and GitHub CLI; follow the repository instructions for release `v1.0.0` assets and checksum verification. From the prepared checkout, use the documented commands:

   ```bash
   python3 scripts/prepare_inputs.py --check-only
   python3 scripts/verify_dataset.py
   python3 scripts/materialize_tasks.py --tasks WFA-OPS-001 WFA-LAW-001 --output ../WorkflowArena-harbor
   ```

   Replace the illustrative task selection with the requested IDs. If inputs are missing, follow `scripts/prepare_inputs.py` and the release download instructions; media tasks additionally require the documented image-loading step. Do not treat successful assembly as successful execution.

2. **Establish a runnable harness contract.** Require a user-provided compatible runner or a subsequently released official framework; do not invent a WorkflowArena launch command. The [paper, Appendix D](https://arxiv.org/html/2611.12856v1#A4) specifies Codex 0.164.0 and Claude Code 2.1.280, a shared 14,400-second execution deadline, and at most 16 concurrent workers. Native delegation is disabled; workers use the same model, harness, and effort through the recorded Workflow interface. Verify these properties in the actual adapter and document deviations. A generic Harbor runner alone does not demonstrate equivalent orchestration support.

3. **Protect the evaluation boundary.** Follow the [evaluation documentation](https://github.com/Einsia/WorkflowArena/blob/main/evaluation/README.md): build solver images from `environment/`, and never mount the full repository or evaluators into a running solver. For tasks using a shared verifier environment, terminate the solver before exposing grading files. Keep judge credentials in the verifier. Readiness fails if evaluator isolation or required grading dependencies cannot be demonstrated.

4. **Capture auditable execution.** Preserve workflow programs, assignments, worker IDs, timestamps, tool events, revisions, terminal artifacts, and verifier logs. Count Workflow starts once per call ID; exclude status/list/read calls. Workflow calls per task includes zero-delegation tasks; worker invocations per Workflow divides total worker invocations by launched Workflows. Use an undefined value when that denominator is zero. These definitions follow [paper Section 4.4](https://arxiv.org/html/2611.12856v1#S4.SS4). Validate trace completeness before computing either statistic, and label any additional diagnostic quantities as local analyses.

5. **Score the integrated deliverable.** Continuous score is `100 * reward`; full resolution requires the domain's complete acceptance criteria, evaluated before rounding. Compute domain means first, then equally average the eight domains for Overall Score and Overall Resolved. Use the released acceptance rules rather than a universal score threshold. An incomplete evaluator is a grading error, not a valid failed submission. These rules are documented in the [evaluation README](https://github.com/Einsia/WorkflowArena/blob/main/evaluation/README.md). Report subset coverage explicitly; never present a selected-task average as the official full-benchmark Overall result.

6. **Separate failure causes and compare reproducibly.** Link each diagnosis to a trace event, artifact, or grading message. Useful local categories include missing worker context, conflicting edits, dependency ordering, unintegrated components, deadline exhaustion, provider errors, and verifier faults; these are diagnostic labels, not additional benchmark metrics. Match task IDs, inputs, runner version, resources, effort, and stopping rules across comparisons. Preserve failed runs and declared exclusions instead of selectively rerunning weak results. Report repeat counts and uncertainty when repetitions exist; do not infer reliable rankings from small differences in single trials.

7. **Control orchestration-transfer experiments.** Hold the executor fixed while changing the workflow source. In static mode freeze the program before execution; in dynamic mode permit revision from worker feedback. The source cannot edit task files, and generation consumes the execution budget. Follow [paper Section 4.5](https://arxiv.org/html/2611.12856v1#S4.SS5). Treat any altered subset or protocol as a separate experiment. Compare against the same executor's self-authored dynamic workflow, and avoid attributing a changed harness or budget to orchestration quality.

## Inputs / Outputs

**Inputs**

- Requested mode: readiness, existing-run analysis, or full execution.
- Benchmark commit, selected task IDs, verified assets, and container environment.
- Compatible runner and adapter, model identifiers, effort, credentials, and execution budget.
- For analysis: complete traces, artifacts, evaluator output, and pricing snapshot if estimating cost.

**Outputs**

- A readiness record within 15 minutes: `ready`, `blocked`, or `analysis-complete`, with checked prerequisites, missing dependencies, and the next executable action. This is an operational status, never a benchmark score.
- A reproducibility manifest recording task/input hashes, commit, image digests, runner and harness versions, configuration, exact commands, and deviations.
- When execution or existing records permit: task-level scores and resolution flags, domain summaries, delegation statistics, resource totals, and evidence-linked failure analysis.

For full execution, preserve the documented task budget rather than shortening it to meet the skill's readiness target. If the runner, assets, credentials, or trace records are unavailable, return the concrete blocker and leave unsupported metrics unset. Aggregate main-agent and worker usage; retain provider token categories and pricing dates so cost calculations can be audited. Do not substitute published aggregate values for measurements from a new run.

## References

- [WorkflowArena: Benchmarking Executable Workflow Orchestration in LLM-Based Multi-Agent Systems, arXiv:2611.12856v1](http://arxiv.org/abs/2611.12856v1).
- [Full paper: protocol, metrics, and execution configuration](https://arxiv.org/html/2611.12856v1).
- [Official WorkflowArena repository and preparation instructions](https://github.com/Einsia/WorkflowArena).
- [Released evaluator documentation and acceptance rules](https://github.com/Einsia/WorkflowArena/blob/main/evaluation/README.md).
- [Official project overview](https://lab.einsia.ai/workflowarena/home.html).
