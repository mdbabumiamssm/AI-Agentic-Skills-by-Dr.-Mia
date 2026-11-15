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
name: 'waza-agent-skill-evaluation'
description: 'Evaluate agent workflows, custom agents, and skills with Waza using representative scenarios, deterministic checks, model-based graders, inspectable evidence, and CI regression gates.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Waza Agent and Skill Evaluation

## Overview

Use Waza to turn agent and skill requirements into repeatable evaluation tasks with explicit graders and retained evidence. Separate deterministic validation from model-based judgments so regressions can be investigated and CI decisions have a clear basis. An offline mock run validates evaluation wiring; only an appropriately configured live run supplies evidence about the evaluated agent.

## When to Use This Skill

- Establish an evaluation suite for a skill, custom agent, or repository workflow.
- Reproduce a reported failure and add a regression case.
- Compare a candidate change against an accepted baseline under the same evaluation conditions.
- Verify output artifacts, tool behavior, and semantic quality with suitable graders.
- Add CI thresholds and preserve evidence for failed as well as successful runs.

## Core Capabilities

1. **Define representative cases.** Map each important requirement to an observable outcome and a stable task ID. Include normal success, boundary inputs, missing prerequisites, malformed inputs, and cases where the correct behavior is clarification or refusal. A negative scenario passes when the agent handles it correctly; do not confuse expected rejection with a failed evaluation.

2. **Select and verify the execution path.** Commands below are grounded in the v0.41.0 documentation linked in References. Record the installed version and consult its help before using flags; main-branch documentation can describe unreleased features. Keep the target revision, fixture revision, executor, model, trial count, and retry policy in the run record.

   ```bash
   waza --version
   waza run --help
   waza gate --help
   ```

   Use an installed, pinned standalone release appropriate to the environment. For a new evaluation of an existing skill in a Waza workspace, `waza new eval my-skill` creates a scaffold; for a repository scenario without a SKILL.md, use `waza new eval inventory --scenario --template repository`. Inspect the generated schema and target binding before editing. For custom agents, follow the custom-agent guide rather than assuming that skill-name binding applies.

3. **Build a bounded suite.** Define the suite in `eval.yaml`, cases under `tasks/`, and input context under `fixtures/`, following the pinned example schema. Use `tasks` globs and `config.trials_per_task`, `config.timeout_seconds`, `config.max_attempts`, and `config.executor` deliberately. Keep evaluator reference answers outside agent-visible fixtures. Freeze the rubric and acceptance policy before measuring the candidate; retain held-out cases when tuning against failures.

   For a first smoke run, use a small representative subset, one trial, and `max_attempts: 1`. Expand trial counts for variability assessment; retries are recovery attempts, not independent successful samples. Budget setup, execution, and inspection within 15 minutes, including an outer job deadline. If the full suite cannot fit, report the completed subset and remaining coverage without claiming full qualification.

4. **Choose graders by evidence type.** Prefer deterministic checks for mechanically decidable requirements and model judgments for semantic criteria that cannot be expressed reliably as assertions.

   | Requirement | Grader choices | Interpretation |
   | --- | --- | --- |
   | Required or forbidden strings and patterns | `text` | Deterministic; matching words alone does not establish correctness. |
   | Structured output, files, or workspace differences | `json_schema`, `file`, `diff` | Deterministic checks over available output and captured artifacts. |
   | Executable assertions or external validators | `code`, `program` | Deterministic only if the validator and dependencies are deterministic. |
   | Recorded tool use and behavioral limits | `tool_calls`, `tool_constraint`, `behavior` | Checks observed events; missing telemetry cannot prove an action never occurred. |
   | Semantic correctness or instruction adherence | `prompt` | Model-based judgment; retain rubric, judge model, and reasons. |

   Example deterministic grader fragment for a task that must return a report marker:

   ```yaml
   graders:
     - type: text
       name: report_marker
       config:
         contains_cs:
           - "REPORT_READY"
         regex_not_match:
           - "(?i)fatal error"
   ```

   This fragment tests formatting only. Add artifact or semantic checks for actual task completion. Waza's `contains` is case-insensitive; use `contains_cs` when case matters. For a `prompt` grader, put the rubric in `config.prompt`, choose an available judge model, and instruct it to call `set_waza_grade_pass` or `set_waza_grade_fail` with reasons. No grade tool call is a failure; a nonzero partial score does not necessarily mean the grader passed. Challenge each grader with known-good and intentionally defective outputs before trusting its verdict.

5. **Preflight and execute with explicit evidence.** With `eval.yaml` and `fixtures/` prepared, run:

   ```bash
   waza preflight eval.yaml --context-dir fixtures --format json > preflight.json
   waza run eval.yaml --context-dir fixtures --executor mock --output smoke-results.json
   ```

   Preflight does not establish credentials, live readiness, or agent quality; ordinary warnings can still exit successfully. A mock pass checks the harness only. When live evaluation is in scope and authentication and model access are available, select the intended executor and a concrete available model:

   ```bash
   waza run eval.yaml --context-dir fixtures --executor copilot-sdk \
     --model "$WAZA_EVAL_MODEL" --output results.json --snapshot snapshots/
   ```

   Set `WAZA_EVAL_MODEL` to the selected model before running; live execution and prompt grading can incur provider usage. Retain results, snapshots, available transcripts and file evidence, grader reasons, command exit status, and version/configuration metadata. Use fresh output locations to prevent stale artifacts from satisfying CI. Fixture copies protect original fixtures but do not provide host filesystem or network sandboxing.

   Inspect every failed or incomplete trial and confirm that required tasks and graders actually ran. Distinguish execution errors, timeout, missing evidence, deterministic assertion failures, and judge disagreement. `waza grade eval.yaml --results results.json` grades saved results without rerunning the task, but prompt graders may still invoke a model; regrading is not a new agent trial.

6. **Enforce regression policy in CI.** Compare an accepted baseline with the candidate using the same suite, fixture inputs, model configuration, graders, and trial/retry policy. Deliberate changes to those conditions require a separately justified comparison or baseline update. Never replace a baseline merely because a candidate failed.

   ```bash
   waza gate --baseline baseline.json --current results.json \
     --max-regression-pct 0 \
     --on-new-tasks fail --on-removed-tasks fail \
     --format github-actions
   ```

   Here zero is an example strict policy, not a benchmark claim: `--max-regression-pct` measures the allowed drop in success rate in percentage points. Select a different tolerance only from the project's stated acceptance policy. Mark critical task definitions `golden: true`; the default golden policy requires them to pass. Default gate exits are 0 for pass, 1 for regression or task-set policy failure, 2 for golden failure, and 3 for configuration error. These differ from `run` and from explicitly selected release-policy modes.

   Preserve producer failures as well as gate failures. Require fresh, parseable, complete results and enforce any absolute pass-rate or per-criterion minimum separately; a relative gate alone can accept an inadequate baseline. Upload available evidence even on failure, then fail the CI job for missing evidence, unsuccessful execution, or a rejected gate. Report per-task changes alongside aggregate rates so improvements elsewhere do not conceal a critical regression. Ordinary gate results do not establish statistical confidence.

## Inputs / Outputs

**Inputs:** Target skill, custom agent, or workflow and its revision; representative prompts and fixtures; observable success and failure criteria; Waza version and execution configuration; grader definitions; trial and retry budgets; and, for regression assessment, an accepted baseline and explicit threshold policy. Live execution additionally requires configured authentication and model access.

**Outputs:** A versioned evaluation specification and task set; a requirement-to-case mapping; preflight diagnostics; fresh run JSON and available snapshots or transcripts; deterministic verdicts separated from model judgments; and a CI gate report with its exit status. Include a concise assessment listing coverage, failures, unavailable evidence, and reproduction commands. A blocked live run or an offline-only smoke result must be labeled accordingly.

**Completion check:** Within the intended 15-minute smoke workflow, validate the suite, execute the declared subset, inspect its evidence, and report a reproducible result. For CI qualification, also verify baseline comparability and enforce the chosen gate; do not count incomplete or mock-only evidence as live agent qualification.

## References

- [Waza source repository](https://github.com/microsoft/waza).
- [Waza v0.41.0 release](https://github.com/microsoft/waza/releases/tag/v0.41.0).
- [Pinned evaluation example and schema usage](https://github.com/microsoft/waza/blob/v0.41.0/examples/code-explainer/eval.yaml).
- [Pinned grader reference](https://github.com/microsoft/waza/blob/v0.41.0/docs/GRADERS.md), [text grader](https://github.com/microsoft/waza/blob/v0.41.0/docs/graders/text.md), and [prompt grader](https://github.com/microsoft/waza/blob/v0.41.0/docs/graders/prompt.md).
- [Pinned CLI reference and gate semantics](https://github.com/microsoft/waza/blob/v0.41.0/site/src/content/docs/reference/cli.mdx).
- [Pinned CI integration guide](https://github.com/microsoft/waza/blob/v0.41.0/site/src/content/docs/guides/ci-cd.mdx).
- [Custom-agent evaluation guide](https://microsoft.github.io/waza/guides/custom-agents/).
