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
name: 'deer-flow-agent-orchestration'
description: 'Set up and operate DeerFlow for long-horizon research, coding, and artifact creation using configured models, sandboxes, skills, memory, and scoped subagents.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# DeerFlow Agent Orchestration

## Overview

Use DeerFlow to carry research, coding, and creation tasks through planning, tool execution, and delivery of inspectable artifacts. This skill covers reproducible setup, provider configuration, scoped delegation, and recovery using the implemented harness rather than assuming that a final chat message proves completion.

Supported baseline: DeerFlow 2.0 source snapshot `c277c1572e5e4c330328a3813adb44601d80d79d`, verified against the repository documentation. This is a documentation pin, not a claim of runtime certification; the legacy `main-1.x` framework is outside scope.

## When to Use This Skill

- Bootstrap or diagnose a DeerFlow installation.
- Configure a model provider, sandbox, reusable skill, or optional tool integration.
- Execute a research, coding, or document-production task with explicit deliverables.
- Delegate independent work while retaining a lead agent responsible for verification.
- Recover an interrupted conversation and inspect existing outputs before retrying.
- Connect a supported messaging channel when the user specifically requests that integration.

## Core Capabilities

1. **Pin and prepare the runtime.** Use the documented 2.0 snapshot and preserve existing installations.

   For a new checkout, use:

   ```bash
   git clone https://github.com/bytedance/deer-flow.git
   cd deer-flow
   git checkout --detach c277c1572e5e4c330328a3813adb44601d80d79d
   git rev-parse HEAD
   ```

   For an existing checkout, inspect its revision and working-tree changes first. Do not force a checkout or overwrite configuration. If another revision is required, read its corresponding README, installation guide, and configuration template before using these instructions.

   Confirm the root contains `Makefile`, `backend/`, `frontend/`, and `config.example.yaml`. Use the lockfiles shipped with the selected revision. Docker requires a reachable daemon and Compose v2.24 or later. Local development requires Python 3.12+, Node.js 22+, pnpm, uv, and nginx; let `make check` identify missing prerequisites.

2. **Configure providers and execution boundaries.** Select an available model with tool-use support, and configure only capabilities needed by the task.

   For a fresh interactive setup, `make setup` launches the provider and execution-preference wizard. For agent-driven or manual setup, use `make config` only when `config.yaml` is absent, then edit the relevant entries from the pinned template. Ensure `models` includes the chosen provider adapter, model identifier, and any required endpoint and environment-variable references. Keep credential values out of reports and do not open secret-bearing files merely to verify them.

   Use `make doctor` to diagnose configuration. A missing credential is a blocker to model execution, not evidence that setup succeeded. Report the required variable name and smallest corrective action.

   Prefer an AIO container sandbox for tasks executing generated code. The local sandbox maps files to host directories and is not a secure execution boundary; host bash is disabled by default. Do not enable host execution merely to bypass a sandbox failure.

3. **Install and launch at the requested scope.** Distinguish preparation from a running application.

   When Docker is available:

   ```bash
   docker info
   docker compose version
   make docker-init
   ```

   When the request includes starting or validating the application, continue with:

   ```bash
   make docker-start
   make docker-logs
   ```

   Without Docker, run `make check`, then `make install` if prerequisites pass; launch with `make dev` when runtime execution is requested. The documented interface is `http://localhost:2026`. Do not describe `make docker-init` alone as successful service startup. Stop on a failed installation step and report its actionable error rather than repeatedly changing dependencies.

4. **Execute a bounded task and verify delivery.** Record the objective, input files, expected artifact, acceptance criteria, and time or cost budget before starting.

   On an already configured deployment, use this smoke task:

   > Create a Markdown checklist explaining the supplied task requirements. Save it as /mnt/user-data/outputs/checklist.md and present the file. Use only the supplied information; do not browse or delegate.

   Run it through the web interface, or use the embedded client from the installed backend Python environment:

   ```python
   from deerflow.client import DeerFlowClient

   client = DeerFlowClient()
   response = client.chat(
       "Create a Markdown checklist with three items: setup verified, "
       "artifact inspected, limitations reported. Save it to "
       "/mnt/user-data/outputs/checklist.md and present the file. "
       "Do not browse or delegate.",
       thread_id="deerflow-smoke-001",
   )
   print(response)
   ```

   Choose a fresh thread ID for an independent test; deliberately reuse the same ID for continuation. Canonical IDs contain 1–64 ASCII letters, digits, hyphens, or underscores.

   Inspect the resulting file through the artifact interface or documented file-access mechanism. Success requires a completed run, a readable nonempty artifact containing the requested items, and a report of the thread identity and artifact location. A returned response object alone does not satisfy these checks.

   Apply the 15-minute outcome to this bounded verification workflow on a prepared environment. Dependency downloads and arbitrary long-horizon tasks may take longer; report incomplete work honestly when the budget expires.

5. **Delegate scoped work and manage reusable context.** Delegate only when parallel latency, specialist capability, or context isolation justifies it.

   Give each child a self-contained task, relevant evidence, output location, and stopping condition. Keep dependent tasks sequential and avoid concurrent writes to the same files. Ordinary DeerFlow `task` calls support `context_mode="isolated"` by default and `"snapshot"` for retained parent conversation context; snapshot mode does not reconstruct messages already removed by compaction. Durable `batch_task` items still need self-contained prompts.

   Inspect child results and underlying artifacts before synthesizing the final answer. Use Capability Center to select skills and plugins relevant to the task. Treat memory as supporting context, not proof that a previous operation succeeded; queued memory extraction can be lost on process exit.

6. **Handle artifacts, messaging, and interrupted runs.** Preserve task identity and existing work before attempting recovery.

   DeerFlow exposes `/mnt/user-data/uploads/` for inputs, `workspace/` for working files, and `outputs/` for deliverables. These are runtime paths, not necessarily host paths. Use `present_files` within the agent workflow to expose final artifacts, then verify they remain accessible after reopening the conversation.

   If the browser or stream disconnects, inspect the existing run before submitting another message. A lost connection does not prove the run stopped. Rejoin or inspect it through the UI or pinned API documentation; an existing-run GET stream is observation-only. Stream-route cancellation requires POST and the appropriate `runs:cancel` permission.

   After an interrupted task, inspect saved files and completed tool actions, then continue in the same thread with an explicit account of what remains. Do not replay external side effects blindly or promise automatic resume after a process restart. Recovery depends on the configured persistence backend and surviving runtime data.

   Use logs and `make doctor` to distinguish provider errors, missing dependencies, sandbox failures, and application failures. Stop retries when the same unresolved cause recurs, preserving the thread/run identifiers and completed outputs.

   Configure messaging only from the pinned README's IM Channels section and matching configuration guide. Require an explicit destination and authorization before sending external messages; installing or configuring a gateway does not itself authorize delivery.

## Inputs / Outputs

**Inputs:** A task objective and acceptance criteria; a new or existing DeerFlow checkout; a provider/model selection with credentials supplied through the configured mechanism; optional input files, skills, tools, memory preferences, and messaging destination; execution and budget constraints.

**Outputs:** The actual repository revision and setup level reached; configuration or dependency blockers without secret values; thread/run identifiers when available; inspected artifact paths or download links; verification results; and an explicit complete, incomplete, or blocked status with the next corrective action.

Do not claim installation, delegated execution, artifact creation, or recovery was tested unless the corresponding operation was performed and checked.

## References

- [Source finding: bytedance/deer-flow](https://github.com/bytedance/deer-flow).
- [Pinned README: runtime, delegation, artifacts, memory, and embedded client](https://github.com/bytedance/deer-flow/blob/c277c1572e5e4c330328a3813adb44601d80d79d/README.md).
- [Pinned installation instructions](https://github.com/bytedance/deer-flow/blob/c277c1572e5e4c330328a3813adb44601d80d79d/Install.md).
- [Pinned configuration template](https://github.com/bytedance/deer-flow/blob/c277c1572e5e4c330328a3813adb44601d80d79d/config.example.yaml).
- [Configuration reference for provider and integration details](https://github.com/bytedance/deer-flow/blob/c277c1572e5e4c330328a3813adb44601d80d79d/backend/docs/CONFIGURATION.md).
- [API reference for thread/run lifecycle and recovery operations](https://github.com/bytedance/deer-flow/blob/c277c1572e5e4c330328a3813adb44601d80d79d/backend/docs/API.md).

