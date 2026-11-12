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
name: 'deer-flow-long-horizon-agent-harness'
description: 'Set up and operate DeerFlow for extended research, coding, and creation tasks using configured models, tools, sandboxes, memory, and controlled subagent delegation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# DeerFlow Long-Horizon Agent Harness

## Overview

Operate DeerFlow as a harness for research, coding, and content creation that require multiple execution steps. This skill connects environment setup, model and tool selection, sandbox boundaries, and artifact verification so extended tasks produce inspectable results. It targets the current harness generation; check the selected revision because the legacy 1.x framework uses a different architecture.

## When to Use This Skill

- Bootstrap or diagnose a DeerFlow environment before starting an extended task.
- Configure models, search tools, skills, and execution permissions for a specific assignment.
- Run research or coding work that benefits from a persistent workspace and bounded delegation.
- Inspect generated files, task failures, or memory behavior before continuing a session.

## Core Capabilities

1. **Revision-aware environment setup.** Read the selected checkout's README, installation guidance, and configuration examples before changing it. Record its revision and reuse an existing installation when possible. The documented quick start uses `make setup`; Docker startup uses `make docker-init` followed by `make docker-start`. The local alternative uses `make check`, `make install`, and `make dev`. Choose one route and verify its prerequisites rather than mixing startup modes. Consult the [repository quick start](https://github.com/bytedance/deer-flow).

2. **Model and tool configuration.** Use the setup wizard or the revision's configuration template to select an available provider, exact model ID, and optional search provider. Keep credentials in the supported environment configuration and exclude their values from reports. Enable only tools and skills needed for the assignment, including explicit bash and file-write preferences. Run `make doctor` to diagnose configuration problems; resolve authentication or tool-call failures before launching extended work. See the [configuration guide](https://github.com/bytedance/deer-flow/blob/main/docs/guide/configuration.md).

3. **Sandbox boundary selection.** Verify the actual execution provider, accessible paths, mounts, and network access before submitting code. Local execution runs on the host; Docker and Kubernetes-backed execution offer different isolation arrangements. Docker socket access can give the gateway host root-equivalent powers, so a container label alone does not establish the intended boundary. Stage only approved inputs and keep generated files within the configured workspace. Consult the [sandbox guide](https://github.com/bytedance/deer-flow/blob/main/docs/guide/sandbox.md).

4. **Bounded task execution.** Define the requested artifact, acceptance criteria, allowed inputs, time budget, and stopping condition before submitting the task through an available documented interface. Start with a small end-to-end task, inspect its result, then expand to the longer assignment. The 15-minute outcome is a target for this bounded verification workflow on a prepared environment, not a promise that installation or hour-scale tasks finish within that period. Report missing prerequisites or incomplete work honestly.

5. **Controlled delegation.** Use subagents only when authorized and useful for independent subtasks. Check that delegation is enabled; the documented web interface associates it with Ultra mode. Set supported runtime caps, such as `subagents.max_total_per_run`, according to the installed configuration schema. Give each child a deliverable, relevant evidence, and distinct write responsibility. Isolated context is the documented default; snapshot context adds retained history but does not restore already-compacted messages. Treat child outputs as material to verify before integration. See the [subagent guide](https://github.com/bytedance/deer-flow/blob/main/docs/guide/subagents.md).

6. **Memory and continuity inspection.** Inspect the configured memory backend and user/agent scope before relying on recalled information. Use Settings > Memory where available to inspect persisted facts; distinguish recalled preferences from current task evidence. Memory extraction can lag behind a completed turn, so do not treat an immediate missing update as proof of data loss. Preserve explicit task notes and artifact locations for handoff. See the [memory guide](https://github.com/bytedance/deer-flow/blob/main/docs/guide/memory.md).

7. **Artifact validation and recovery.** Open the actual generated files and compare them with acceptance criteria. For research, check cited sources; for code, run relevant checks within the approved environment. Record the session or run identifier when available, artifact paths, validation results, and unresolved errors. On failure, distinguish provider, tool, sandbox, and task-content problems; correct the identified cause before a bounded retry. Preserve partial work and avoid rerunning completed external actions blindly.

## Inputs / Outputs

**Inputs**

- A DeerFlow checkout or running instance and its selected revision.
- A concrete task, source files or URLs, expected output format, and acceptance criteria.
- Available model credentials supplied through supported configuration, required tools, and optional skills.
- Approved workspace paths, execution mode, network permissions, memory policy, delegation allowance, and resource budget.

**Outputs**

- A concise environment record: revision, model identifier, enabled tools, execution mode, and diagnostic status without secrets.
- Verified task artifacts with accessible paths and evidence of the checks performed.
- A completion or partial-progress report identifying remaining work, failure causes, and the next actionable step.

**Bounded acceptance example:** On a prepared instance, provide a short local text file and request a summary saved in the permitted output directory. Confirm that the resulting file exists, opens correctly, reflects the input, and respects the agreed tool and delegation boundaries. If any condition fails or the budget expires, report the observed state instead of declaring success.

## References

- [DeerFlow repository and source finding](https://github.com/bytedance/deer-flow).
- [Model, tool, and environment configuration](https://github.com/bytedance/deer-flow/blob/main/docs/guide/configuration.md).
- [Sandbox execution and filesystem guidance](https://github.com/bytedance/deer-flow/blob/main/docs/guide/sandbox.md).
- [Subagents and delegation controls](https://github.com/bytedance/deer-flow/blob/main/docs/guide/subagents.md).
- [Long-term memory behavior](https://github.com/bytedance/deer-flow/blob/main/docs/guide/memory.md).

These links track upstream `main`; use the corresponding files at the selected checkout revision when behavior differs.
