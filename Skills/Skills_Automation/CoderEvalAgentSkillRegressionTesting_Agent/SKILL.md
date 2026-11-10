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
name: 'coder-eval-agent-skill-regression-testing'
description: 'Test skills, MCP servers, and CLIs with coding agents using Coder Eval sandboxed YAML tasks, observable success criteria, controlled A/B experiments, and CI release gates.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Coder Eval Agent Skill Regression Testing

## Overview

Use Coder Eval to test whether a coding agent can discover and use a skill, MCP server, or CLI to complete representative tasks. Define observable outcomes in sandboxed YAML suites, compare baseline and candidate configurations, and use the evidence to identify regressions before a release. This workflow evaluates agent behavior and produced artifacts; it does not establish general model rankings or evaluate multi-agent orchestration design.

## When to Use This Skill

- A skill's instructions, examples, dependencies, or discovery metadata change.
- An MCP server or CLI changes its interface, installation, or output format.
- A model, agent runtime, prompt, or tool configuration needs comparative evaluation.
- A release needs repeatable CI checks for task completion and skill activation.
- A refresh pipeline needs to investigate failures after environment or model drift.

## Core Capabilities

1. **Establish an executable evaluation contract.** Identify the target revision, intended agent, representative user requests, baseline, and release requirements. Check the installed `coder-eval --version`, agent runtime, credentials availability, and matching upstream documentation before choosing YAML fields or flags. Record prerequisites that are missing instead of treating an unrun evaluation as a pass.

2. **Construct representative tasks.** Select a normal workflow, a meaningful boundary case, and a previously observed failure when available. For skill discovery, use natural user requests that do not explicitly name the skill; keep explicit-invocation tests separate. Include a negative routing case when accidental activation matters. Use small, deterministic fixtures with expected outputs, and avoid exposing grading answers in the agent prompt.

3. **Define observable success criteria.** Author task YAML using the installed version's task schema, including `task_id`, `initial_prompt`, agent configuration, sandbox settings, and `success_criteria`. Combine artifact checks with executable behavioral assertions; a file existing or a command exiting successfully alone may not prove correctness. Use `skill_triggered` where supported to assess activation independently of output quality. Set scoring weights, thresholds, and mandatory release conditions before running; prevent partial credit from masking a critical failure.

4. **Control environment and tool access.** Pin the framework, agent runtime, target revision, dependencies, fixtures, and sandbox image where feasible. Record the exact model identifier, tool permissions, prompt configuration, resource limits, and credential source without recording secrets. Configure skill installation and MCP/CLI availability according to the selected harness documentation, then verify they are visible inside the sandbox. A temporary directory is not a security boundary: use the documented container driver for untrusted execution and keep production data and credentials out of fixtures.

5. **Validate and run a bounded smoke evaluation.** Start with a small task that exercises the actual integration, then inspect both its artifact and criterion results. With the CLI already installed and the task authored, use the following sequence, substituting the actual task and recorded run paths:

   ```bash
   coder-eval plan tests/tasks/skill_smoke.yaml
   coder-eval run tests/tasks/skill_smoke.yaml
   coder-eval report runs/latest
   ```

   Planning checks task configuration without spending model tokens; it does not prove runtime success. Configure supported execution limits before running. Treat 15 minutes as the smoke-workflow budget, not a promise that a complete regression suite will finish; report incomplete coverage if the budget expires. Preserve the specific run directory because `runs/latest` can change after another run.

6. **Compare baseline and candidate fairly.** Follow the version-matched A/B experiment schema and vary only the intended factor. For a skill revision comparison, hold tasks, graders, agent settings, dependencies, and environment constant while changing the skill revision. Run variants in fresh workspaces, record repetitions, and retain every attempt. Compare per-task outcomes, activation evidence, elapsed time, and available token/cost telemetry; report missing measurements explicitly. A single stochastic result is provisional evidence, not a reliable improvement claim.

7. **Diagnose failures and set CI gates.** Separate invalid task definitions, setup/authentication failures, missing tool exposure, activation failures, incorrect outputs, grader defects, and timeouts. Inspect traces and artifacts before changing instructions or thresholds. Apply the predefined release policy to mandatory checks and comparative results, and retain failing evidence for reproduction. For GitHub Actions, use the documented action inputs for the pinned revision, provision the agent runtime separately, preserve its failure exit status, and archive JUnit/run reports even on failure. Keep credentials unavailable to untrusted fork execution; redact reports before publishing them.

## Inputs / Outputs

**Inputs**

- Target skill, MCP server, or CLI with baseline and candidate revisions.
- Representative prompts, fixture data, expected behavior, and any existing failure reproductions.
- Agent/runtime/model settings, sandbox requirements, dependency pins, and authorized test credentials.
- Evaluation budget, repetition policy, mandatory checks, and acceptable regression thresholds.

**Outputs**

- Validated YAML tasks and, when requested, an A/B experiment configuration or CI gate configuration.
- A reproducibility record containing versions, resolved settings, fixture revisions, commands, and exact run paths.
- Per-task criterion results and supporting artifacts/traces, with activation distinguished from task success.
- A baseline/candidate comparison containing attempt counts, observed failures, available resource measurements, and limitations.
- A release recommendation of pass, fail, or inconclusive under the predefined policy, with focused remediation for failures.

**Completion check:** The smoke task has actually run, its expected artifact has been inspected, its criteria have been evaluated, and a report identifies the target revision and outcome. Configuration-only validation, unavailable credentials, interrupted execution, or missing required evidence must be reported as incomplete. A broader release recommendation requires the agreed regression coverage; a smoke pass alone is insufficient.

## References

- [UiPath/coder_eval](https://github.com/UiPath/coder_eval) — source finding and framework overview.
- [Task Definition Guide](https://github.com/UiPath/coder_eval/blob/main/docs/TASK_DEFINITION_GUIDE.md) — YAML schema, criteria, and scoring.
- [A/B Experiments](https://github.com/UiPath/coder_eval/blob/main/docs/AB_EXPERIMENTS.md) — comparative configuration workflow.
- [Run-Limit Parity](https://github.com/UiPath/coder_eval/blob/main/docs/agents/HARNESS_PARITY.md) — harness-specific interpretation of execution limits.
- [Docker Isolation](https://github.com/UiPath/coder_eval/blob/main/docs/DOCKER_ISOLATION.md) — container execution configuration.
- [CI Gate and GitHub Action](https://github.com/UiPath/coder_eval/blob/main/docs/CI_GATE.md) — CI integration and reporting.
