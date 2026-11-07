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
name: 'ml-research-artifact-reproducibility'
description: 'Reconstruct ML experiments from papers, code, and supporting artifacts, preserve procedural decisions, and audit execution evidence for reproducible handoffs.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# ML Research Artifact Reproducibility

## Overview

Turn scattered experiment instructions into a traceable reproduction plan, tested execution path, and evidence-based assessment. ReproBridge provides the methodological reference: preserve a refined codebase, environment setup, experiment launchers, and a reproduction guide so later attempts can reuse recovered knowledge. The operational checklist below is a local adaptation, not the released ReproBridge implementation.

## When to Use This Skill

- A paper's results depend on undocumented configuration or setup choices.
- Released code, supplementary material, and manuscript instructions disagree.
- A successful experiment needs a reusable handoff for another researcher.
- Reported reproduction scores need verification against execution evidence.
- A compute-limited attempt needs an honest readiness or blocker report.

## Core Capabilities

1. **Verify source availability.** Read the complete method and relevant appendices; locate official artifacts through primary links. The full ReproBridge paper was checked on 2026-11-06; its project page still advertised code and benchmark releases as forthcoming. Recheck availability before adopting implementation-specific commands, package names, or schemas; never guess a repository URL.

2. **Define a bounded task.** Assign a task ID and identify the exact table cell or claim, dataset version, split, model, metric units, and training versus evaluation scope. Record the available hardware, access restrictions, and time budget. Group shared computation without merging distinct reported measurements.

3. **Maintain a decision ledger.** Use one entry per ambiguity with these fields: `decision_id`, `task_id`, `source_locator`, `evidence_class`, `gap`, `proposed_resolution`, `rationale`, `validation`, and `status`. Evidence classes are `documented`, `observed_in_code`, `agent_inferred`, and `unresolved`. A successful run does not convert an inference into a documented fact. Keep explicit contradictions separate from missing details; for paper reproduction, prefer explicit paper requirements and record conflicting code.

4. **Prepare a controlled execution.** Inspect setup scripts before running them. Pin the source revision and record local diffs, dependency versions, accelerator/runtime versions, and input checksums where feasible. Use an isolated working environment within the authorized workspace. Check file availability, split membership, and configuration parsing before expensive computation. Record repairs as execution fixes, scientific corrections, or convenience edits, with affected tasks and checks.

5. **Capture a measurable attempt.** Fix seeds, repetitions, aggregation, and acceptance criteria before examining results. Assign a unique run ID; retain the command, working directory, configuration snapshot, start/end times, exit status, raw output locations, and metric-extraction procedure. Keep diagnostic runs separate from scored runs. Never replace measurements with bundled scores or repeat runs solely to obtain a favorable result.

6. **Audit the evidence chain.** Trace each reported value back through extraction, raw outputs, executed code, and task specification. Check split leakage, incomplete training, wrong metric units, stale output reuse, and unsupported substitutions. Record each finding with a file/line or output locator, its possible measurement effect, and counterevidence. Separate execution completion, numerical agreement, and integrity judgments; a mismatch alone is not an integrity violation. Preserve uncertainty when evidence is insufficient.

7. **Test the handoff.** Provide the four artifact components with explicit coverage and limitations. Attempt a fresh-shell rerun using only the handoff; record whether full execution, a small real experiment, or merely a smoke test was completed. Independent review and fresh-session reproduction strengthen the assessment; if unavailable, label the result self-reviewed and independent validation pending. Historical outputs remain evidence, never substitutes for new execution.

## Inputs / Outputs

**Inputs:** Paper/version and selected claims; official repository revision; configuration and supplementary materials; data/checkpoint locations; existing logs; resource limits; authorized output location; and any predefined comparison rule. Missing inputs become explicit blockers rather than inferred downloads or invented values.

**Outputs:** A task manifest, decision ledger, environment record, scoped patch record, runnable entry points, reproduction guide, run evidence index, and task-level audit report. These are logical deliverables: combine them into one report if the user permits only a single file, and create executable artifacts only within the requested scope.

Use this compact result schema:

```text
task_id | run_id | execution_status | observed_value | metric_units
reference_source | comparison_rule | numerical_status
evidence_locators | integrity_status | unresolved_decisions | next_action
```

Use execution statuses `not_run`, `blocked`, `smoke_only`, `partial`, or `completed`; numerical statuses `not_assessed`, `compatible`, or `mismatch`; and integrity statuses `not_reviewed`, `supported`, `uncertain`, or `invalid`. Explain every non-completed task. Without a justified comparison rule, report the observed difference and leave numerical compatibility unassessed.

**15-minute completion contract:** Spend approximately three minutes on intake, four on specification and gaps, five on a feasible check, and three on evidence review and handoff. Deliver a valid scoped assessment within that window; full ML training may require a separately budgeted continuation. Valid output requires source locators, decisions, coverage, evidence or blocker reasons, and a next action for each selected task. Do not claim scientific reproduction from setup success, a shortened workload, or an unexecuted plan.

## References

- [ReproBridge: Building Reproducibility into ML Research Artifacts, arXiv:2611.08452v1](http://arxiv.org/abs/2611.08452v1) — supplied source finding.
- [Full method and appendices](https://arxiv.org/html/2611.08452v1) — Section 3 and Appendix B describe artifact construction; Section 4.1 and Appendix C describe execution auditing.
- [ReproBridge project page](https://wyyyb.github.io/reprobridge_website/) — release availability checked on 2026-11-06; code and benchmark links were listed as forthcoming.
