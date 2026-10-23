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
name: 'adaptive-skill-granularity-routing'
description: 'Route coding tasks among primitive tools, macro skills, and hybrid execution, then recover from observable macro failures while preserving workspace state and verifying results.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Adaptive Skill Granularity and Recovery

## Overview

Choose an execution granularity for a coding task and revise execution after observable macro failures. This skill makes fallback inspectable and bounded while retaining useful work and requiring evidence that the requested behavior is correct.

## When to Use This Skill

- A coding task can be performed through either individual tool actions or a reusable composite skill.
- A retrieved macro appears relevant, but its applicability to the current workspace is uncertain.
- A macro reports failure after potentially modifying files.
- An executor needs explicit fallback rules and an auditable execution trace.
- A team is evaluating whether macro savings justify failure and recovery costs.

## Core Capabilities

1. **Establish the evidence boundary.** The cited paper uses structural pre-execution features, an offline depth-limited decision tree, and a finite-state executor with one-way macro-to-micro downgrade. Its reported experiments use small deterministic task sets and controlled faults; successful late-failure states can make unconditional recovery wasteful. These findings support testing adaptive execution locally, not assuming production-wide gains. The operational rules below are practical adaptations, not a reproduction of its learned policy.

2. **Define completion before choosing tools.** Identify the requested behavior, affected workspace, existing edits, acceptance checks, and time or action budget. Record whether checks can actually run. Treat the 15-minute outcome as a bounded workflow target, not a guarantee: if verification remains unavailable or unsuccessful at the limit, return an incomplete result with the remaining work.

3. **Assess each macro as an executable contract.** Read its instructions and record its arguments, applicability conditions, permitted changes, completion checks, observable failure signals, and fallback entry point. Confirm that required tools exist and that the macro stays within the user's scope. A descriptive match alone does not establish compatibility. If side effects cannot be determined, inspect further or choose primitive actions.

4. **Select and record the initial route.** Use evidence available before execution, including file types, dependency structure, operation type, and macro compatibility. Do not use the current run's eventual verifier result as a routing input. The following decision rules are local defaults, not thresholds claimed by the paper:

   | Route | Choose when | Execution behavior |
   | --- | --- | --- |
   | Micro | No macro meets the contract, or the change requires exploratory inspection | Inspect, edit, and check through individual actions. |
   | Macro | The operation fits a known contract and stopping on failure is acceptable | Execute the macro; stop on failure unless a previously defined recovery policy permits continuation. |
   | Hybrid | The macro fits and an inspectable, bounded repair path exists | Start with the macro; permit primitive recovery after a recorded failure. |

   Record the selected route, recovery permission, rationale, macro identity, and remaining budget. Prefer micro execution when compatibility evidence is insufficient. Use a learned router only when its model and feature schema are available and locally evaluated; otherwise disclose that selection is heuristic.

5. **Recognize failures through evidence.** Define the trigger before invoking the macro. Useful local triggers include a missing replacement target, nonzero tool exit, rejected arguments, invalid output structure, failed syntax or acceptance checks, and an executor exception. Record the actual observation and whether changes occurred. Distinguish a precondition failure with no edits from a failure that leaves partial changes; do not infer either state solely from an error message.

6. **Recover from the observed workspace.** Use this local transition procedure:

   - `READY → MICRO_ACTIVE` or `READY → MACRO_ACTIVE` follows the recorded route.
   - A completed action enters `VERIFY`; only passing the task's required checks permits `DONE`.
   - A macro failure enters `ASSESS`: inspect affected files, compare against the recorded baseline, and run applicable checks before making further edits.
   - If the requested final state already passes all required checks, enter `DONE` and record a macro error without crediting a recovery.
   - Otherwise, if recovery was enabled and state is understood, log `ASSESS → MICRO_ACTIVE` before issuing the next action. Repair only the unmet requirements and verify again.
   - If recovery is disabled, state cannot be safely reconstructed, or the budget is exhausted, enter `STOPPED` with the failure evidence and remaining requirements.

   Preserve valid partial changes and pre-existing user edits. Avoid replaying a module move or replacement that already occurred. Undo only identified task-owned changes when necessary for repair; never reset the whole workspace as an automatic fallback. After downgrade, keep the run in micro execution to avoid repeated macro failure loops. Stop a retry sequence when it repeats the same failure without new evidence or a changed repair action.

7. **Verify both output and control flow.** Check the final workspace against the original acceptance criteria, not a weaker substitute chosen after failure. Inspect the change scope and report check failures or unavailable checks. For recovery, verify that an observed trigger preceded the downgrade and that the continuation used the preserved state. A useful scenario review covers a compatible macro, an inapplicable macro, a partial edit requiring repair, an error after completed work, and a failure with recovery disabled. Run such cases in an authorized disposable environment when evaluating an executor.

8. **Measure without overstating success.** Record actions, elapsed time, observed token or monetary costs when available, macro attempts, failure events, transitions, and final verification status. Define metric denominators explicitly: macro adoption can mean runs attempting a macro divided by all runs; recovery rate can mean verified successes after actual repair divided by failed runs that needed repair and had recovery enabled. Report already-complete error cases separately. Compare routes under equivalent tasks, starting states, and checks; include failed-run costs and mark cost per success undefined when no run succeeds. Route agreement measures conformity to a specified policy, not correctness by itself.

## Inputs / Outputs

**Inputs**

- Task requirements, workspace location, current file state, and applicable repository instructions.
- Available primitive tools and candidate macro instructions or contracts.
- Acceptance checks, execution budget, and any existing failure observations.
- Optional historical run data and a validated router artifact with its feature schema.

**Outputs**

- Initial route and rationale, including whether fallback is enabled.
- An execution record containing actions, observed failures, state assessments, and transitions.
- Requested code changes, preserved valid partial work, and final verification evidence.
- A terminal status of verified success or incomplete/stopped, with unresolved requirements and measured costs where available.

Keep records in the response unless the task authorizes persistent artifacts. This skill does not itself authorize additional files, external actions, or expansion of the coding task.

## References

- Dingxuan Zhang, Yanhui Jiang, Jun Wang, and Hefeng Zhou. [Adaptive Skill Granularity for Coding Agents: Runtime Routing and Macro-to-Micro Recovery](http://arxiv.org/abs/2610.31461v1), 2026. [Full text](https://arxiv.org/html/2610.31461v1), Sections 4–5: routing and recovery mechanism, evaluation design, and limits of unconditional fallback. The full text was checked when authoring this skill; no experimental results were independently reproduced.
