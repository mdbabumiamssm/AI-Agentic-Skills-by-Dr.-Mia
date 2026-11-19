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
name: 'langfuse-agent-observability-and-evaluation'
description: 'Instrument OpenTelemetry-based agent traces with Langfuse, build evaluation datasets, attach scores, and compare prompt or harness revisions to diagnose failures and regressions.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Langfuse Agent Observability and Evaluation

## Overview

Use Langfuse to connect agent execution traces with evaluation evidence and versioned prompts. This skill produces a reproducible comparison of agent behavior, helping locate tool failures, quality regressions, and latency changes on Langfuse Cloud or an existing self-hosted instance. Langfuse provides OpenTelemetry-based observability alongside evaluation and prompt management; this workflow focuses on those capabilities ([repository](https://github.com/langfuse/langfuse)).

## When to Use This Skill

- Instrument an agent, retrieval step, tool call, or model generation for debugging.
- Convert representative, sanitized requests into repeatable evaluation cases.
- Attach deterministic checks, human ratings, or rubric-based judge scores to executions.
- Compare prompt versions, model settings, or agent harness revisions before release.
- Diagnose missing traces, broken parent-child relationships, or incompatible SDK integrations.

## Core Capabilities

1. **Select compatible SDKs.** Inspect the application language, dependency lockfile, deployment endpoint, and server version before changing instrumentation. The official documentation checked on 2026-11-19 lists Python SDK v5 and JS/TS SDK v6 for Langfuse server v4.59.0 or later; self-hosted v3 deployments need an older compatible SDK. Recheck the compatibility matrix at execution time and record exact installed versions. Do not mix legacy trace APIs with current observation APIs ([SDK overview](https://langfuse.com/docs/observability/sdk/overview)).

2. **Capture an interpretable execution tree.** Create a root observation for each request, with nested observations for agent steps, retrieval, tools, and generations. Preserve OpenTelemetry context across asynchronous work; prefer the application's existing integration when it already records the same calls. Record sanitized inputs and outputs, errors, model identity, and available usage data. Confirm a representative trace reaches the intended project before expanding coverage.

3. **Construct evaluation datasets.** Define cases with stable identifiers, inputs, optional expected outputs, and metadata describing the scenario. Include ordinary requests and relevant failure cases; do not treat generated expectations as verified ground truth. Keep reference answers out of the agent's task inputs. Python provides `create_dataset` and `create_dataset_item`; preserve the dataset version or an immutable local snapshot for repeatability ([datasets](https://langfuse.com/docs/evaluation/experiments/datasets)).

4. **Attach interpretable scores.** Specify each metric's name, type, range or categories, direction, and evaluator revision before execution. Prefer deterministic checks for exact requirements; use documented rubrics and sampled human review for subjective judgments. With Python, `create_score` accepts `trace_id`; observation-level scores also require the corresponding `observation_id`. The current low-level `api.scores.create` accepts a list, so do not substitute an old single-score payload ([scores via SDK](https://langfuse.com/docs/evaluation/evaluation-methods/scores-via-sdk)).

5. **Compare revisions reproducibly.** Use `run_experiment` with a task callback and local data, or the documented hosted-dataset runner. Current SDKs removed older dataset-run helpers; verify signatures before adapting examples. Run baseline and candidate against the same cases and evaluators, retaining per-item results and failures ([experiments via SDK](https://langfuse.com/docs/evaluation/experiments/experiments-via-sdk)). Pin prompt versions for comparisons and record the resolved version if selecting by label, because labels can move ([prompt management](https://langfuse.com/docs/prompt-management/get-started)).

### Execution workflow

1. **Establish the smallest useful run.** Identify the agent entry point, baseline and candidate revisions, evaluation question, and available time and model-call budget. Use a small representative subset for the 15-minute smoke workflow; a smoke result is not a production quality benchmark. If credentials or a runnable harness are missing, report the blocker and provide an unexecuted integration plan without claiming success.

2. **Configure the existing project.** Read `LANGFUSE_PUBLIC_KEY`, `LANGFUSE_SECRET_KEY`, and `LANGFUSE_BASE_URL` from the existing environment or secret store. Use the user's selected region or self-hosted endpoint. Never include credentials in traces or reports. Sanitize sensitive content before export, including tool arguments and exception details. Keep external tools in a test environment or use fixtures when evaluating actions that would alter real systems.

3. **Instrument and verify one request.** For current Python SDKs, initialize with `from langfuse import get_client`, then `langfuse = get_client()`. Use `langfuse.start_as_current_observation(as_type="span", name="agent-request")` as a context manager and nested observations for subordinate work; use `as_type="generation"` for actual model calls. In JS/TS, initialize OpenTelemetry and its Langfuse span processor before instrumented application code. Flush the SDK or span processor before a short-lived process exits. These entry points are documented in the [SDK overview](https://langfuse.com/docs/observability/sdk/overview).

4. **Freeze cases and evaluation rules.** Store the dataset identity, case identifiers, expected-output provenance, evaluator configuration, and scoring policy. Distinguish task failures, evaluator failures, and missing observations. Failed or unscored cases must remain visible in the final denominator rather than disappearing from the comparison.

5. **Execute both revisions.** Keep model settings, tool fixtures, and sampling policy constant unless they are the intended change. Record the harness commit, resolved prompt version, model configuration, and run identifier. Start with bounded concurrency; respect the configured call budget. For transient failures, retry only within a finite budget and preserve attempt metadata so retries do not become independent successes.

6. **Confirm ingestion and compare.** Check that representative traces, nested observations, and scores are retrievable in the intended project. Successful local execution or a flush call alone does not prove ingestion. Compare paired cases, score coverage, task errors, and observed latency; compare token usage or cost only where recorded and comparable. Inspect regressions individually and report uncertainty for small or stochastic samples. Leave missing telemetry as missing rather than zero.

7. **Return evidence and a decision.** Summarize baseline versus candidate against user-defined acceptance criteria, linking the experiment and representative traces. If no release threshold exists, describe the differences without declaring a release pass. Prompt promotion, production rollout, and infrastructure deployment are separate from running this evaluation.

### Troubleshooting and completion checks

- **No telemetry:** check credentials and region, compatible server/SDK versions, exporter configuration, application initialization order, and process shutdown. Inspect SDK/exporter errors without exposing secrets; do not repeatedly replay costly requests merely to test connectivity.
- **Disconnected or duplicated spans:** check context propagation and overlapping automatic/manual instrumentation. Maintain one logical execution tree per evaluated request.
- **Missing scores:** verify trace and observation identifiers, score types, and ingestion status. Use stable score identifiers where supported when retrying delivery to avoid duplicate measurements.
- **Invalid comparisons:** rerun only affected cases when versions, case selection, or evaluator settings differ; otherwise label the results as non-comparable.
- **Completion:** provide compatible version details, dataset identity, two distinguishable revision runs when comparison was requested, retrievable trace and score evidence, and a concise result with failures and limitations. If the time budget expires first, return partial evidence and the exact remaining blocker.

## Inputs / Outputs

**Inputs**

- Agent source or callable harness; language, dependencies, and deployment context.
- Existing Langfuse project endpoint and credentials supplied through secure configuration.
- Sanitized cases or permission to derive them from specified traces; expected outputs or evaluation rubrics.
- Baseline and candidate prompt/harness revisions, model settings, and optional acceptance criteria.
- Execution budget, concurrency limits, and permitted test-tool environment.

**Outputs**

- Scoped instrumentation changes or a clearly labeled integration plan if execution is blocked.
- Dataset identifier/version or local snapshot identity, with case schema and expectation provenance.
- Experiment identifiers and trace links, score definitions, and evaluator configuration.
- Comparison report containing case counts, score coverage, errors, measured performance where available, regressions, and limitations.
- Explicit execution status: verified, partially verified, or blocked; never invented trace IDs, scores, or benchmark claims.

## References

- [Langfuse source repository](https://github.com/langfuse/langfuse) — supplied source finding and platform scope.
- [Langfuse SDK overview and compatibility](https://langfuse.com/docs/observability/sdk/overview) — current instrumentation entry points and server requirements.
- [Langfuse datasets](https://langfuse.com/docs/evaluation/experiments/datasets) — dataset and case creation.
- [Langfuse experiments via SDK](https://langfuse.com/docs/evaluation/experiments/experiments-via-sdk) — task runners, evaluation, and migration notes.
- [Langfuse scores via SDK](https://langfuse.com/docs/evaluation/evaluation-methods/scores-via-sdk) — score types and trace/observation association.
- [Langfuse prompt management](https://langfuse.com/docs/prompt-management/get-started) — versioned prompts and labels.
