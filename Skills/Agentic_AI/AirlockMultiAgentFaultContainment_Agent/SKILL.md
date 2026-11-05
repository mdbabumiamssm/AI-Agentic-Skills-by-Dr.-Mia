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
name: 'airlock-multi-agent-fault-containment'
description: 'Reproduce and evaluate selective post-alert containment of multi-agent faults using narrowed permissions, dependency provenance, quarantine, and guarded cuts while measuring healthy work retention.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Airlock Multi-Agent Fault Containment

## Overview

Use this skill to evaluate whether an alert stops faulty outputs from causing downstream external actions while unaffected work continues. Airlock combines EdgeCitadel's parent-to-child permission narrowing with an append-only provenance DAG and two containment tiers: full quarantine and an optimized cut. The workflow below is a bounded reproduction procedure, not a claim to reproduce the paper's complete evaluation in fifteen minutes.

## When to Use This Skill

- A faulty artifact has already reached consumers in other agent branches.
- Source cancellation leaves downstream tasks able to act on contaminated inputs.
- A workflow needs selective containment without restarting all healthy tasks.
- An evaluation must distinguish post-alert leaks, convergence, and work retained without redo.
- A proposed recovery operation needs evidence that releasing a hold preserves containment.

## Core Capabilities

1. **Verify the implementation before execution.** Read the versioned paper, especially Sections 3–4 and Appendix A.6, and inspect the selected repository revision. Record actual controller entry points, switch APIs, provenance schema, test commands, and recovery support. The main README describes runtime installation; it alone does not establish availability of the Airlock evaluation harness. Do not invent commands or substitute a current package installation for a pinned research artifact.

2. **Audit authority boundaries.** Check each child scope against both its parent scope and administrator grant. Attempt a forbidden target and excessive delegation in a disposable fixture; require observable denial. Inspect where credentials and external effects are mediated rather than trusting an agent's declared permissions.

3. **Track cross-branch influence.** Connect the faulty output version to its actual readers, subsequent shared-memory turns, pending actions, and later consumers. Keep the permission tree distinct from the provenance DAG. Compare the computed affected set with independently specified fixture ground truth.

4. **Exercise both containment tiers.** Full quarantine freezes affected work without requiring a solver. Evaluate the optimized cut separately for its ability to reduce healthy-work disruption; never infer correctness from a lower cut cost alone.

5. **Verify release and recovery boundaries.** The paper distinguishes reversible task, egress, and event holds from permanent grant, admission, and credential revocations. Its gateway serializes effect checks and controls under one lock. Optimized convergence requires fenced confirmations, journal ingestion through those fences, source-side labels, reconciled taint, and separation from protected sinks. Release only temporary holds whose removal preserves separation; revocations stay permanent. Missing evidence is indeterminate. Airlock does not repair memory or undo committed effects. These constraints come from Sections 3–4 and Appendix A.6; use them as acceptance conditions, not assumed implementation features.

6. **Measure containment and useful progress independently.** Report faulty commits and healthy work completed without rerun from destination observations and task identities. Continued unprotected output from a tainted agent is not proof that its answers became correct.

### Reproduction workflow

Use an already prepared sandbox and mock external sinks by default. Execute only commands confirmed in the selected revision. The following time allocations are a smoke-test budget, not measured runtime guarantees.

**Minutes 0–3: establish evidence and scope.**

- Record paper version, repository commit, runtime configuration, and dependency availability.
- Map the verified implementation to the six capabilities above. Mark absent mechanisms explicitly.
- Identify protected sink classes, fault origins, alert timing, and the planned test endpoint.
- If the artifact or necessary instrumentation is unavailable, produce a blocked report with the missing evidence; do not present a standalone model as an Airlock reproduction.

**Minutes 3–6: prepare a minimal deterministic fixture.**

- Have a producer emit a corrupted version that a consumer in another branch reads before the alert.
- Make that consumer schedule a protected mock action and propagate the input to a later task or shared-memory turn.
- Run an independent healthy branch that never consumes the corrupted version, with its own observable output.
- Record task IDs, artifact versions, consumption events, alert timestamp, scheduled deadlines, control acknowledgments, and destination commits.
- Define ground-truth affected tasks independently of the controller's own taint labels.

**Minutes 6–11: run paired containment checks.**

- Reuse the same seed, plan, injected fault, and alert point for source-only cancellation, full quarantine, and optimized containment when supported.
- Include an early alert and an alert after consumption. Allow the fixture to reach its declared endpoint so that delayed leaks cannot disappear through premature observation cutoff.
- Exercise control delay, late journal arrival, and a newly spawned consumer after the alert. Attempt a protected effect during the control transition.
- Inspect permission denials and whether healthy task IDs finish without restart. Preserve the independent branch as a positive progress control.

**Minutes 11–14: inspect convergence and guarded release.**

- For the optimized tier, locate the actual convergence evidence and timestamp; a command acknowledgment alone is insufficient.
- Attempt a release that would reopen a tainted path and require refusal. Verify an eligible temporary hold can be removed while the protected path remains blocked.
- Keep holds when confirmation or journal evidence is missing. Report unresolved convergence rather than treating a timeout as successful containment.
- Check persistent shared state explicitly; starting a new invocation must not silently make contaminated memory clean.
- Treat memory repair, replay, new credentials, or compensating physical actions as separate recovery work requiring its own verified procedure.

**Minutes 14–15: return the evaluation record.**

For each policy and scenario, report:

- Protected faulty commits from alert to endpoint, and separately from verified convergence to endpoint.
- Convergence latency when established; otherwise an explicit unknown value and reason.
- Healthy work retained without redo: healthy task weight completed without rerun divided by eligible healthy task weight. State weights and denominator; use unavailable when the denominator is zero.
- Healthy tasks stopped, affected tasks contained, rejected permission attempts, and release decisions with evidence pointers.

A clean post-alert result requires zero faulty protected commits throughout the entire post-alert window and completion of the observation period. A zero post-convergence count does not erase earlier leaks. Fail the safety check for a protected leak after verified convergence or an unsafe release; mark incomplete evidence indeterminate. Claim healthy-work preservation only for the tasks actually observed, and distinguish this smoke test from the publication's broader results.

## Inputs / Outputs

**Inputs:** Versioned paper and pinned repository or prepared checkout; documented harness commands; delegation scopes and grants; task/artifact identities; journal and dependency records; fault/alert fixture; protected sink policy; timing assumptions; and independent destination observations.

**Outputs:** A concise report containing the source-to-implementation mapping, executed commands, fixture configuration, permission results, affected-set comparison, containment actions, convergence evidence, release checks, and policy comparison metrics. Include remaining revocations, unresolved recovery needs, and a pass, fail, blocked, or indeterminate outcome with reasons. Keep raw observations available within the user's authorized output scope.

**Limits:** The paper's guarantee assumes mediated effects, trusted enforcement components, and a correct initial fault-origin set. Missing provenance, bypassing the gateway, or unresolved timing evidence prevents an assurance claim; a fifteen-minute report can validly document these limitations without claiming a successful reproduction.

## References

- [Airlock: Bounding the Multi-Agent Blast Radius via Selective Online Containment, arXiv:2611.06727v1](http://arxiv.org/abs/2611.06727v1) — source finding and reported evaluation.
- [Airlock full text](https://arxiv.org/html/2611.06727v1) — runtime design, convergence and release conditions, recovery limits, and evaluation definitions.
- [EdgeCitadel repository](https://github.com/EdgeCitadelTeam/EdgeCitadel) — inspect and pin the actual implementation before executing experiments.
- [EdgeCitadel README](https://github.com/EdgeCitadelTeam/EdgeCitadel/blob/main/README.md) — runtime setup and operational documentation; mutable reference, not a pinned experiment specification.
