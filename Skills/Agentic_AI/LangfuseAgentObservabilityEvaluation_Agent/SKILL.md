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
name: 'langfuse-agent-observability-evaluation'
description: 'Instrument agents with Langfuse OpenTelemetry tracing, inspect execution failures, run dataset evaluations, and compare prompt versions with reproducible evidence.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Langfuse Agent Observability and Evaluation

## Overview

Use this skill to connect agent execution traces to repeatable quality evaluations and prompt-version comparisons in Langfuse. Langfuse supports OpenTelemetry-based observability, evaluations, and prompt management on Cloud or self-hosted deployments; this workflow turns those capabilities into inspectable evidence for debugging and regression decisions.

## When to Use This Skill

- Add tracing to an agent, model call, retrieval step, or tool invocation.
- Diagnose missing observations, incorrect parent relationships, slow steps, or execution errors.
- Evaluate agent outputs against a fixed dataset and explicit scoring criteria.
- Compare two prompt versions using the same inputs and model configuration.
- Adapt an integration to the installed Langfuse SDK and server versions.

## Core Capabilities

1. **Select compatible integration APIs.** Read the application manifest and lockfile, identify the runtime, framework, SDK version, and Cloud or self-hosted endpoint, then consult the matching official SDK documentation. At authoring verification on 2026-11-27, the documentation specifies Python SDK v5 and JS/TS SDK v6 for Langfuse v4, with a self-hosted minimum of v4.59.0; older servers need the documented compatible SDK. Recheck compatibility before changing dependencies. Use the existing framework integration when available, and avoid mixing legacy trace APIs with current observation APIs. See [SDK overview and compatibility](https://langfuse.com/docs/observability/sdk/overview).

2. **Instrument meaningful execution boundaries.** Create an active root observation for one agent request, with child observations for retrieval, tool execution, and model generations. Preserve context across asynchronous work and avoid duplicate instrumentation of the same model call. Current Python APIs include `get_client()`, `observe()`, and `start_as_current_observation()`; JS/TS uses tracing helpers such as `startActiveObservation` with an initialized OpenTelemetry exporter. Record errors, model identity, and available usage without inventing missing token or cost data. See [instrumentation guidance](https://langfuse.com/docs/observability/sdk/instrumentation).

3. **Inspect trace delivery and structure.** Execute a synthetic request, finish observations, and flush the SDK or shut down the OpenTelemetry exporter before a short-lived process exits. Locate the trace in the intended project and verify parent-child links, final output, error status, and timings. If it is missing, check credentials, endpoint, exporter initialization, sampling, context propagation, and flush behavior before repeating execution. Redact credentials and sensitive payloads before export; local execution does not imply local-only telemetry.

4. **Run dataset-based experiments.** Define a task that consumes each item's input and returns the actual application output. Use a small fixed dataset with stable item identifiers, expected outputs where appropriate, and explicit evaluator semantics. Current Python supports `langfuse.run_experiment(...)`; JS/TS supports `langfuse.experiment.run(...)`. Confirm callback signatures against the installed SDK; removed dataset-run helpers must not be copied into newer integrations. Configure bounded concurrency, timeouts, and a run name. Hosted and local datasets are supported, but local experiment data is still exported to the configured instance. See [experiments via SDK](https://langfuse.com/docs/evaluation/experiments/experiments-via-sdk).

5. **Score outputs and compare prompt versions.** Prefer deterministic checks for schema validity, exact answers, and tool constraints; use rubric-based judges or human review for qualitative criteria. Record score name, type, range, direction, evaluator version, and failure handling. Keep evaluator failures separate from low task scores. Retrieve explicit prompt versions, compile their variables, and attach prompt identity to the relevant generation using the version-matched API. Evaluate baseline and candidate on identical dataset items with fixed model settings; report per-item regressions as well as aggregate scores. See [prompt management](https://langfuse.com/docs/prompt-management/get-started).

### Bounded execution workflow

- **Minutes 0–3: establish inputs.** Select the requested tracing, evaluation, or comparison task. Confirm project endpoint, installed versions, available credentials, a runnable agent entry point, and a small synthetic or approved dataset. Read credentials from the environment, never from committed examples. Use `LANGFUSE_PUBLIC_KEY`, `LANGFUSE_SECRET_KEY`, and the SDK's documented base-URL setting; current SDK documentation uses `LANGFUSE_BASE_URL`.
- **Minutes 3–7: verify one trace.** Instrument the smallest useful path, run one request, flush, and inspect the exported trace. Exercise an error case if error tracing is part of the task. Do not count a successful local function return as evidence of telemetry delivery.
- **Minutes 7–12: measure the requested behavior.** Run the bounded dataset experiment or two prompt variants. Record dataset identity or content hash, code revision, prompt versions, model settings, evaluator configuration, and run identifiers. Use sandboxed tools or recorded fixtures when evaluation would otherwise repeat real-world side effects.
- **Minutes 12–15: report evidence.** Provide trace and experiment links or identifiers, completed/failed item counts, scoring denominators, observed regressions, and limitations. Distinguish measured latency and usage from unavailable metrics. If prerequisites or service access prevent completion, return the exact blocker and partial evidence; do not claim a successful run.

Success means the requested path has verifiable evidence: a visible trace with correct nesting for instrumentation, persisted item results and scores for evaluation, or comparable baseline and candidate results for prompt comparison. The 15-minute target assumes an accessible project and runnable application; it is a workflow budget, not a platform performance guarantee. A small smoke dataset verifies integration and does not establish general model quality.

## Inputs / Outputs

**Inputs**

- Application entry point, runtime/framework, dependency versions, and code revision.
- Langfuse project endpoint and environment-provided credentials; provider credentials when real model calls are required.
- Synthetic or approved sample requests, dataset identity, expected outputs, and evaluator rubric or scoring functions.
- Prompt names and explicit baseline/candidate versions for comparisons.
- Execution budget, concurrency limit, payload-redaction requirements, and task-specific acceptance thresholds when supplied.

**Outputs**

- Minimal integration changes or a concrete configuration proposal when execution is blocked.
- Trace identifiers or links with verified observation structure and delivery status.
- Experiment identifiers, per-item outputs/scores/errors, and aggregate results with denominators.
- A reproducibility record covering versions, dataset, model settings, prompts, evaluators, and code revision.
- A concise findings report identifying observed regressions, missing evidence, and the next corrective action; no unsupported accuracy claims or automatic production prompt promotion.

## References

- [Langfuse source repository](https://github.com/langfuse/langfuse) — supplied source finding for the platform and its scope.
- [Langfuse SDK overview](https://langfuse.com/docs/observability/sdk/overview) — runtime setup, SDK/server compatibility, and export lifecycle.
- [Langfuse instrumentation](https://langfuse.com/docs/observability/sdk/instrumentation) — active observations, decorators, and trace context.
- [Langfuse experiments via SDK](https://langfuse.com/docs/evaluation/experiments/experiments-via-sdk) — datasets, task callbacks, evaluators, and experiment runners.
- [Langfuse prompt management quickstart](https://langfuse.com/docs/prompt-management/get-started) — prompt retrieval, compilation, and versions.
