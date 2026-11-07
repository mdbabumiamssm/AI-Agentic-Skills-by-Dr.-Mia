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
name: 'executable-research-environments'
description: 'Turn research papers into executable environments with reproduced baselines, open research objectives, fair evaluation, and traceable artifacts using the Scientist World workflow.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Executable Research Environments for AI Scientists

## Overview

Convert a research paper into a runnable setting where new methods can be compared with a reproduced reference method. Follow the reproduce, design, package, and verify sequence described in Scientist World, keeping the research question separate from its proposed solution. Use the released implementation for concrete interfaces and preserve evidence for every reproduction or improvement claim.

## When to Use This Skill

- A paper needs an executable baseline before exploring extensions.
- An AI scientist needs a bounded research objective, fixed evaluation rules, and usable data.
- A released Scientist World environment needs preparation, inspection, or validation.
- Research results need an auditable record connecting methods, runs, scores, and source artifacts.

## Core Capabilities

1. **Select a feasible research question.** Identify the central empirical claim, accessible data, executable method, and resource requirements. Prefer a question that admits multiple solutions and measurable progress; defer candidates whose essential experiments cannot run with available resources.
2. **Resolve and pin released artifacts.** Follow links from the versioned paper to its code and environment. Record repository commits, package versions, data revisions, licenses, and retrieval dates; inspect actual files before choosing commands.
3. **Reproduce the reference method.** Run the paper's method and relevant simple baselines under the same evaluation conditions. Separate published numbers, locally reproduced numbers, and results from reduced experiments.
4. **Define an open evaluation contract.** Specify inputs, outputs, metric direction, fixed resources, allowable changes, and visible versus held-out evaluation. Keep scoring independent of method implementation.
5. **Package repeatable execution.** Capture dependencies, data acquisition, configurations, seeds, and evaluator entry points so a fresh workspace can reproduce the measured result.
6. **Verify and preserve provenance.** Recompute scores from saved outputs, audit leakage and constraints, and retain successful and failed experiments. Mark proposed extensions explicitly until executed and evaluated.

### Workflow and acceptance gates

**1. Establish scope and a 15-minute checkpoint.** Record the paper version, chosen claim, available hardware, execution budget, and requested mode. Within 15 minutes, return an environment assessment, evaluation contract, artifact manifest, and either a completed bounded smoke test or a precise blocker with the next executable step. This checkpoint does not certify full reproduction; do not invent results when downloads, training, or verification exceed the time budget. Stop or leave jobs running only according to the agreed execution scope.

**2. Inspect the release before implementation.** The system repository separates paper selection, environment construction, and evaluation; the standalone `sciworld` CLI consumes existing environments. An arbitrary paper ID is not evidence that an environment exists. Find a matching release first; otherwise treat construction as a new task. The published builder stack requires cluster infrastructure, so do not imply that installing the CLI provisions it. See the system and CLI references below.

For an installed, configured CLI, these documented commands support inspection and preparation. Replace `PAPER_ID` with a confirmed environment ID and use the installed version's help to confirm options:

```bash
sw list --search PAPER_ID --json
sw select PAPER_ID
sw info PAPER_ID --full
sw research PAPER_ID --mode improve --agent codex --setup-only
```

`--setup-only` prepares a workspace without launching research. Choose `innovate` when the scientist should receive the problem without the reproduced method. Preserve the untouched source separately from each trial. Container backends provide workspace isolation; the documented `venv` backend does not. The CLI launches unattended agents with permission prompts disabled, so choose an execution boundary consistent with the user's authorized scope.

**3. Write the research contract.** State the question without prescribing the source method. Specify data splits, output schema, metrics and their direction, aggregation, compute limits, and forbidden shortcuts. For an existing environment, retain its fixed contract. For a new one, define tolerances and acceptance criteria before examining candidate improvements. Keep model size, training data, and compute comparable unless changing those quantities is explicitly the research question.

**4. Establish reproduction evidence.** Inspect the paper's implementation and preprocessing, then run a small diagnostic before committing the full budget. Record the exact command, seed, configuration, hardware, exit status, logs, and output hashes. Compare each measured metric with the corresponding published condition; document discrepancies rather than silently changing evaluation. A reduced dataset or compute budget yields an adapted experiment and needs its own reference scores. Confirm that simplification has not removed the scientific challenge.

**5. Build or map the workspace.** The released CLI documents these paths; inspect each environment because its dependencies and experiment commands can differ:

| Artifact | Purpose |
| --- | --- |
| `READMEs/problem.md`, `evaluation.md` | Research question and evaluation interface |
| `READMEs/TARGETS.md`, `CONSTRAINTS.md`, `DIRECTION.md`, `EXPERIMENTS.md` | Objectives, boundaries, permitted changes, and experiments |
| `data/`, `baseline/`, `method/`, `eval/` | Inputs, sanity baseline, candidate or source method, and visible evaluator |
| `scripts/download.sh`, `baseline.sh`, `evaluate_train.sh` | Data preparation and visible execution entry points |
| `environment/container.def` | Software environment recipe |
| `scoring/scores_train.json`, `scores_test.json` | Measured visible and held-out outputs |
| `idea.md`, `reflection.md`, `status.json` | Hypothesis, interpretation, and trial completion state |

For a new environment, use this layout as a compatibility target rather than claiming compatibility from filenames alone. Verify the selected release's parser and output contracts. Pin dependencies and data checksums, document required access, and replace machine-local assumptions with repeatable setup. Keep the hidden evaluator and source construction records outside the scientist's accessible workspace.

**6. Run extensions under frozen rules.** Develop against visible evaluation, preserving each hypothesis, change, and outcome. Freeze the final method before held-out scoring; do not tune against the held-out results. The released launcher introduces `scripts/evaluate_test.sh` in its test phase. Inspect `sw status --json`, then use `sw evaluate PAPER_ID` with the appropriate run selector when several trials exist. A generated score file alone does not establish validity: check evaluator integrity and recompute from outputs in a trusted verification context.

Report raw metrics and signed gains against the reproduced reference: candidate minus reference for higher-is-better metrics, reference minus candidate for lower-is-better metrics. Use official FNI only when the release supplies applicable normalization coefficients; record their revision and date. Do not relabel an ordinary percentage change as FNI or compare FNI across different coefficient sets. Preserve visible and held-out scores separately even if the CLI supplies an aggregate.

**7. Verify before marking ready.** Rerun setup and scoring from a clean workspace, confirm input/output schemas and finite metrics, and check that the method cannot alter fixed evaluation or access held-out material during development. Reconcile the report with logs and score files. If independent verification was not performed, state that explicitly. Failed reproduction, missing data, and incomplete evaluation remain visible blockers rather than passing results.

## Inputs / Outputs

**Inputs:** A versioned paper URL or identifier; source code or implementation details; accessible datasets and weights; a chosen research claim; execution resources and budget; requested `improve` or `innovate` mode; and an output location authorized by the user.

**Outputs:** A prepared environment or a clearly scoped construction specification, a fixed evaluation contract, and a provenance manifest linking the paper, source commit, dependency versions, data checksums, configuration, commands, and run artifacts. Include a results table with metric, split, published value, reproduced value, candidate value, signed gain, and evidence location; use `not run` for absent measurements.

Assign an explicit status: `blocked`, `prepared`, `smoke-tested`, `reproduced`, or `verified`. These are local reporting labels, not Scientist World's official milestone names. Describe the exact reproduction scope and verification performed. Keep proposed extensions separate from executed experiments and reproduced findings.

Retain failures and resource usage in the research record. Before sharing artifacts, exclude credentials, including CLI-created authentication copies and `launch.env`; preserve the licenses of source code, datasets, weights, and paper-derived material. The system's license does not automatically cover every environment asset.

## References

- [Scientist World: Grow a Research World for AI-Scientists, arXiv:2611.08478v1](http://arxiv.org/abs/2611.08478v1) — source finding and conceptual workflow.
- [Scientist World system repository](https://github.com/zlab-princeton/scientist_world) — released construction infrastructure, component boundaries, and requirements.
- [Scientist World CLI documentation](https://github.com/zlab-princeton/scientist_world_cli/blob/main/README.md) — command syntax, workspace layout, execution modes, and scoring behavior inspected for this skill.
- [ScientistWorld environment repositories](https://github.com/ScientistWorld) — discover the environment matching a selected paper; pin the specific repository revision before execution.
