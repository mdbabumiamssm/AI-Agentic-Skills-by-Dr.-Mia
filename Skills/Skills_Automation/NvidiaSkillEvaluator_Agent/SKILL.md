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
name: 'nvidia-skill-evaluator'
description: 'Evaluate AI agent skills with NVIDIA SkillEvaluator using quality gates, semantic overlap checks, synthetic datasets, and live comparisons for skill refresh and regression review.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# NVIDIA SkillEvaluator

## Overview

Use NVIDIA SkillEvaluator to assess skill quality, detect overlapping guidance, generate evaluation cases, and measure behavior with and without a skill. This workflow turns daily skill refreshes into reviewable evidence, separating static findings from observed task outcomes.

## When to Use This Skill

- Review a new or changed skill before adding it to a catalog.
- Investigate duplicated or conflicting instructions across skills.
- Build task cases for a skill whose behavioral impact is untested.
- Compare a refreshed skill with its previous revision and a no-skill baseline.
- Produce a regression report for a scheduled quality-control run.

## Core Capabilities

1. **Reproducible setup.** Record the evaluator revision, skill revision, dependency environment, configuration, and selected stages before running checks.
2. **Quality gates.** Run validation and preserve actionable findings; distinguish failed criteria from unavailable checks.
3. **Semantic overlap detection.** Inspect internal repetition and catalog overlap, retaining the compared catalog revision.
4. **Synthetic evaluation cases.** Generate cases grounded in intended tasks, then review their expected outcomes before using them as evidence.
5. **Live behavioral comparisons.** Compare matched tasks with and without the skill to identify improvements and regressions.
6. **Regression reporting.** Preserve configuration, case-level outcomes, and unresolved findings so later refreshes can be compared fairly.

### Evaluation procedure

**Scope and setup.** Select the skill directory, comparison catalog if relevant, previous report, and acceptance criteria. Choose a workload that fits the available time; the 15-minute target is a workflow budget, not a guarantee that installation or live evaluation will finish. If prerequisites are missing, return an incomplete report identifying the blocker.

Consult the [repository documentation](https://github.com/NVIDIA/SkillEvaluator) for the selected release. Use an isolated checkout pinned to a recorded commit and its lockfile. The documented source-checkout runtime uses Git, uv, and Python 3.13. Confirm command options with `--help` before adapting these examples; replace `./my-skill` with the actual skill directory.

```bash
uv run --python 3.13 skillevaluator quality-check ./my-skill -r cli -r json -o reports/quality
uv run --python 3.13 --extra tier2 skillevaluator tier2 ./my-skill
uv run --python 3.13 --extra tier3 skillevaluator tier3 ./my-skill
```

Run only the requested stages. The quality check can run without a model; overlap analysis needs model services, and live evaluation needs a configured agent and sandbox. A static pass does not establish security coverage or behavioral improvement.

**Interpret findings.** Capture each command, exit code, report location, and finding. Treat missing dependencies, credentials, scanner execution, or timed-out runs as incomplete evidence. For a completed failing gate, report the affected file, criterion, supporting evidence, and proposed correction. Review semantic matches for shared purpose and contradictory advice before recommending consolidation; similarity alone does not justify deleting a skill. Mark catalog comparison as unperformed when no catalog was supplied.

**Prepare cases.** Supply task requests and observable expected outcomes using the release-supported generation interface. Include normal tasks, boundary cases, and requests where the skill should not activate. Review synthetic cases for ambiguity, unsupported assumptions, and answers leaked into prompts. Freeze the reviewed dataset before comparing revisions; record case identifiers and a dataset hash.

**Compare behavior.** Use identical cases, model versions, tool access, scoring rules, and budgets for the no-skill and candidate-skill conditions. Start each condition in a fresh sandbox or session to prevent carryover. For refresh review, also run the previous skill revision when available. Retain task outputs and traces, verify whether the skill was actually loaded, and identify unexpected tool actions. Repeat cases when practical to assess variability; label single-run observations accordingly.

**Report and decide.** Compare candidate and baseline outcomes per case, then summarize aggregate changes with denominators. Keep the no-skill comparison separate from the previous-revision comparison. Apply acceptance criteria chosen before execution and report unmet criteria without relaxing thresholds after seeing results. Preserve partial results at the time limit and list unexecuted stages explicitly. Do not claim improvement when missing runs or changed conditions prevent a valid comparison.

## Inputs / Outputs

**Inputs:** A skill directory containing `SKILL.md` and supporting assets; evaluator revision and configuration; optional catalog and prior skill revision; reviewed task cases and scoring criteria; provider and sandbox configuration for model-backed stages; execution budget and output location. Keep credentials outside reports.

**Outputs:** Available evaluator reports and logs; overlap findings with review decisions; generated and reviewed cases; paired task outcomes and trace locations; and a regression summary containing:

- Evaluator, skill, catalog, and dataset versions or hashes.
- Executed stages, configuration, exit codes, and elapsed time.
- Gate findings and their pass, fail, incomplete, or skipped status.
- Per-case baseline and candidate results, aggregate denominators, and observed variability.
- Acceptance decision, corrective actions, and remaining coverage gaps.

A valid deliverable accurately records the executed scope and evidence. An incomplete evaluation is a reportable outcome, never a passing gate or proof of benefit.

## References

- [NVIDIA/SkillEvaluator — source repository, setup, tier workflows, and reporting documentation](https://github.com/NVIDIA/SkillEvaluator).
