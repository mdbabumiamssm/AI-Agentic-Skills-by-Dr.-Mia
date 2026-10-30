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
name: 'computer-agent-task-verification'
description: 'Verify computer-agent task completion through explicit criteria, adaptive GUI or code checks, and traceable evidence for evaluation, data curation, and training.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Computer Agent Task Verification

## Overview

Verify desktop, mobile, and web task outcomes by gathering evidence for explicit completion criteria instead of accepting an agent's completion claim. Inspired by VeriCUA, this skill uses context inspection, rubric construction, adaptive GUI or code checks, and evidence assessment to address incomplete views of application state.

### Source and implementation status

The versioned full text was inspected on 2026-10-30, including Sections 3 and Appendix A. No official implementation for this exact paper was confirmed in the inspected sources; the similarly named GitHub project in References cites a different paper and must not supply installation instructions or results for this one. This skill is an operational adaptation, not an official software distribution or a claim of reproduced performance.

## When to Use This Skill

- Audit a computer-using agent's claimed task completion.
- Label completed attempts for evaluation or dataset selection.
- Investigate discrepancies between trajectory judgments and application state.
- Assess whether a proposed training reward has sufficient evidence.
- Compare active verification with a fixed, trajectory-only judge.

## Core Capabilities

1. **Context-grounded criteria.** Inspect available context before turning the request into testable requirements.
2. **Adaptive inspection.** Choose GUI observations or code checks that address unresolved criteria.
3. **Evidence capture.** Associate every finding with the actual observation and check that produced it.
4. **Explicit arbitration.** Distinguish confirmed completion, confirmed failure, and insufficient evidence.
5. **Paired evaluation.** Compare outcomes with a trajectory-only baseline on the same attempts.

### Confirmed method rules

The paper distinguishes hard and soft requirements, with met, unmet, or unknown statuses. Hard requirements derive from the request, not discovered storage locations. Findings must cite executed checks. Attested observations can support outcomes; asserted intentions or completion claims cannot establish success. Frame timing and state persistence constrain evidence use.

Online confirmation requires persistent attested evidence; refutation requires two corroborating probes. Offline confirmation requires support from every completed audit path; refutations require separate probes. The default offline configuration allows three paths and early stopping after two agree on every hard requirement. All hard requirements met yields success; any unmet yields failure. Online unresolved tasks abstain. Offline unresolved tasks use an initial judgment withheld from exploration, without changing unknown requirements. See the [method and implementation details](https://arxiv.org/html/2610.40578v1).

### Operational workflow

The timing, report structure, and conservative default below are local workflow choices.

1. **Scope the attempt and access, minutes 0–2.** Record the original request, attempt ID, platform, terminal-state identity, available observations, and permitted inspection interfaces. Determine whether evidence is recorded-only or includes a live/snapshotted environment. Do not assume shell access provides a GUI controller. If an authorized inspection interface is absent, report the access limitation.
2. **Establish the baseline and rubric, minutes 2–4.** If comparison is requested, obtain a trajectory-only verdict before exploration, using a fixed observation bundle and prompt; hide it from subsequent auditing. Inspect context and give each requirement an ID, request excerpt, predicate, and hard/soft designation. Specify whether it concerns an end state or a change; a change requires an appropriate earlier reference.
3. **Investigate, minutes 4–11.** For each unresolved criterion, select the least costly available check likely to answer it. Record the check inputs, returned observation, timestamp or frame index, source identity, errors, and criterion association. Use only authorized inspection operations. Never perform the original task to make verification pass. Avoid writes, submissions, installations, and arbitrary generated shell execution; use approved read operations and constrained predicates. Record any inspection-induced state disturbance.
4. **Arbitrate, minutes 11–13.** Apply the declared evidence policy consistently. A failed tool call is an access error, not proof of task failure. Repeated readings of the same cached artifact do not constitute independent corroboration. Resolve conflicting observations by checking target identity, timing, and scope; retain unresolved conflicts. Default to abstention when evidence remains insufficient. For an explicit paper-reproduction run, enable the documented offline fallback and label it separately as fallback-derived.
5. **Report and compare, minutes 13–15.** Produce the record below, even when blocked or timed out. Stop once the verdict is supported, further inspection cannot resolve the remaining criteria, or the budget expires. Successful execution of this skill means a valid verification report, not necessarily a success verdict for the underlying task.

### Evaluation against trajectory-only judgments

Use identical attempt IDs and independently obtained reference labels. Keep labels hidden during verification; record model, prompt, observation selection, elapsed time, and inspection budget for both conditions. Separate gains from additional state access from gains within the same recording; a live-state comparison does not isolate reasoning quality.

Report sample count, agreement with reference labels, and false positives divided by reference-negative attempts. Also report decision coverage and accuracy among decided attempts, with raw counts and undefined denominators marked unavailable. Preserve abstentions in the report; if a binary score maps them to failure, disclose that mapping. For offline fallback runs, show fallback frequency and evidence-only coverage separately. Without reference labels, report paired disagreements and evidence differences, not accuracy or claimed improvement.

## Inputs / Outputs

**Inputs:** Original task request; agent attempt and completion claim; retained screenshots, action logs, or structured observations; optional authorized terminal-state access and pre-state reference; inspection budget; optional independent evaluation labels. Missing sources must be recorded, not invented.

**Output:** One verification report, returned inline unless artifact storage is authorized, containing:

```yaml
attempt_id: supplied-or-derived-identifier
mode: offline-or-online
policy: conservative-or-paper-reproduction
source_manifest: []
rubric: [] # ID, request basis, hard/soft, predicate, met/unmet/unknown
checks: [] # ID, inputs, source, time, observation/excerpt, error, evidence attributes
findings: [] # requirement ID, check IDs, support/refutation, rationale
verdict: success-or-fail-or-abstain
verdict_basis: evidence-or-offline-fallback
unresolved: []
limitations: []
baseline_comparison: null # Paired results and metrics when available
elapsed_seconds: null
stop_reason: criteria-settled-or-budget-or-access-limit
```

Before returning, ensure every finding references an executed check, every mandatory criterion has a status, the verdict follows the chosen policy, and the report distinguishes observations from inference. Include enough excerpts or artifact locators to review the decision without exposing unrelated sensitive content.

## References

- [VeriCUA: Verification as Agentic Exploration for Computer-Using Agents, arXiv:2610.40578v1](http://arxiv.org/abs/2610.40578v1) — source finding.
- [Versioned full text](https://arxiv.org/html/2610.40578v1) — Sections 3 and Appendix A describe verification and evidence rules.
- [Jeffjeno/VeriCUA](https://github.com/Jeffjeno/VeriCUA) — name-disambiguation reference only; its README identifies a different work, arXiv:2610.26469. It is not a confirmed implementation of the source finding.
