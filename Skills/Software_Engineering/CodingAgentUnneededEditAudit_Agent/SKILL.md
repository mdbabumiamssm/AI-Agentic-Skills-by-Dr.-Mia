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
name: 'coding-agent-unneeded-edit-audit'
description: 'Audit coding-agent patches for test-relative unnecessary edits using trajectory attribution and individual block reversions, separating code necessity from patch size.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Coding Agent Unneeded Edit Audit

## Overview

Audit a resolved repository patch using Delta Residue Attribution (DRA): identify surviving edit blocks, trace their origins, and evaluate necessity by reverting them individually. Larger patches alone do not establish over-editing; report test-relative evidence and retain uncertainty about requirements the tests miss.

## When to Use This Skill

- A reviewer requests evidence that a coding-agent patch contains unnecessary changes.
- A passing patch includes unrelated behavior, refactoring, instrumentation, tests, or documentation.
- A harness evaluation needs code necessity alongside task success.
- A team wants a reviewable pruning proposal with reproducible test evidence.

## Core Capabilities

1. **Establish a reproducible audit boundary.** Record the base revision, submitted patch hash, task requirements, environment, test commands, and available time. Work in an isolated checkout; preserve the submitted patch and any unrelated working-tree changes. Run the full-patch control first. A failing control blocks necessity grading: report the failure rather than interpreting variant results.

2. **Replay and attribute edits.** Replay recorded file mutations in order and require exact agreement with the submitted final state. Attribute a unit to the latest step writing or deleting its lines. Preserve action identifiers and supporting log excerpts. Reconstruct file writes from artifacts without blindly executing historical shell commands. If reconstruction fails, label provenance unknown; continue only as an explicitly limited patch-only audit. Distinguish debug, fix, and cleanup intent using trajectory evidence, with confidence and rationale for ambiguous labels.

3. **Build a stable unit inventory.** DRA merges contiguous change blocks separated by at most one unchanged line; created or deleted files each form one unit. Classify tests, documentation, configuration, nonsemantic edits, and code separately. Assign stable IDs with file paths and base/final ranges. Explain classification ambiguities, including mixed-purpose blocks. Unsupported binary changes or rename handling must be listed as limitations rather than silently omitted. Do not equate an ungraded test or documentation change with a defect.

4. **Grade individual counterfactuals.** From the unchanged full patch, revert exactly one unit and rerun the fixed oracle. Rejection means needed; acceptance means unneeded under that oracle. Restore the full state before every next trial. Freeze test selection and fixtures so removing an edit cannot silently remove its own checks. For SWE-bench, use FAIL_TO_PASS and PASS_TO_PASS when available; disclose when only repository and agent-visible tests are accessible. Test edits reset by the harness remain ungraded.

   Paper replication counts base-only reversions, parse failures, repeated hangs, and evaluation errors as necessary. Record these conventions separately from observed test failures. For practical review, use an explicitly adapted status of inconclusive for infrastructure failures, flaky runs, unappliable variants, or timeouts; do not present that adapted score as exact replication. Allow at most one diagnostic rerun within the budget, then retain uncertainty.

5. **Review candidates against intent.** This operational extension checks each test-unneeded candidate against acceptance criteria and changed interfaces before proposing removal. Preserve evidence for error handling, compatibility, performance, and other requirements absent from the test suite. Optional differential probes should record inputs and outputs for the full and reverted states; a trusted reference patch can provide another comparison. Behavioral difference alone proves neither usefulness nor dispensability. Keep probe conclusions separate from the original test verdict.

6. **Report measured coverage and prune cautiously.** Necessity precision is needed units divided by all units in the stated population; report code separately. For incomplete practical audits, report needed, unneeded, inconclusive, ungraded, and untested counts; label any scored-subset ratio explicitly and use N/A for a zero denominator. Do not claim global minimality. Compare developer patches only with matching units, tests, and environments. If pruning is requested, remove one candidate at a time in isolation, rerun validation after each accepted removal, then validate the combined patch. Independent passing reversions do not guarantee their combination passes.

### Execution budget and completion

Use the first few minutes to verify inputs and time the control run. Reserve time for reporting; estimate remaining trials from measured test duration and prioritize code units relevant to the review question. Stop starting trials when they would exceed the 15-minute budget, terminate overdue processes, and report pending units explicitly. A bounded partial report is valid output but is not a completed full-patch audit.

Completion requires a recorded control result, a unit inventory, reproducible evidence for every assigned verdict, disclosed coverage, and a concrete recommendation. Missing prerequisites produce a blocked report naming the missing artifacts; they do not justify invented results. Deliver findings without changing the user's working patch unless that change is within the requested scope.

## Inputs / Outputs

**Required inputs**

- Repository at a known base revision and the exact final patch.
- Task statement and acceptance criteria.
- Reproducible test environment, commands, fixtures, and pass/fail criteria.

**Inputs for additional capabilities**

- Ordered trajectory with tool arguments, observations, and file snapshots for attribution.
- Developer patch for calibrated comparison; unavailable hidden tests must remain explicitly unavailable.
- Audit scope, runtime budget, and whether a pruning proposal or applied change is requested.

**Outputs**

- Audit status: complete, partial, or blocked, with elapsed time and limitations.
- Manifest identifying revisions, patch, environment, oracle, and replay fidelity.
- Per-unit table: ID, path/ranges, category, writing step/role, reversion status, test result, verdict, evidence location, and recommendation.
- Code and non-code counts, grading coverage, and clearly scoped necessity ratios.
- Candidate removals with requirement-level rationale; optional reduced patch with combined validation results.

Example verdict wording: “U07 is unneeded under the recorded repository tests; compatibility behavior remains unverified, so retain pending review.” Avoid claiming that passing tests prove an edit is unnecessary in every deployment.

## References

- [How Much Do Coding Agents Over-Edit, and Why? An Empirical Study of Unneeded Edits in Repository-Level Agent Patches](http://arxiv.org/abs/2611.34875v1), arXiv:2611.34875v1. Method: Section 3 and Appendix B; behavioral validation: Appendix C. [Full text](https://arxiv.org/html/2611.34875v1). The workflow above adds explicit operational budgets, audit statuses, and review safeguards; these are implementation guidance, not additional empirical findings.
