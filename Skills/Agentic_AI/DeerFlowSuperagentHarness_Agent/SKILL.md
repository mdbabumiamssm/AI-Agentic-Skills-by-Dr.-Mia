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
name: 'deer-flow-superagent-harness'
description: 'Use Deer Flow, ByteDance''s open-source long-horizon SuperAgent harness, to research, code, and create with sandboxes, memories, tools, skills, subagents, and a message gateway.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Deer Flow SuperAgent Harness

## Overview

This skill guides agents in evaluating, configuring, and operating Deer Flow, ByteDance's open-source long-horizon SuperAgent harness for research, coding, and creation tasks. Use it when the work benefits from coordinated sandboxes, memory, tools, skills, subagents, and a message gateway rather than a single short-lived prompt.

## When to Use This Skill

- Building or assessing a long-horizon agent workflow that may run for minutes to hours.
- Comparing Deer Flow with other agentic harnesses, multi-agent frameworks, or deep-research stacks.
- Planning a Deer Flow deployment, local setup, sandbox integration, tool registry, skill layer, memory store, or message gateway.
- Debugging task orchestration across research, coding, artifact creation, and subagent handoffs.
- Turning a broad research or engineering objective into an auditable agent workflow with checkpoints and deliverables.

## Core Capabilities

1. **Repository orientation** - Inspect the Deer Flow source tree, entrypoints, configuration, examples, and documentation before proposing changes.
2. **Harness setup** - Identify required runtimes, environment variables, model/provider settings, sandbox configuration, and startup commands from the repository.
3. **Workflow design** - Break long-horizon tasks into research, planning, execution, review, and delivery stages that match Deer Flow's orchestration model.
4. **Tool and skill integration** - Map user goals to available tools, skill definitions, memory behavior, and subagent roles without inventing unsupported integrations.
5. **Sandbox operations** - Use repository-supported sandbox mechanisms for coding, execution, validation, and artifact generation when available.
6. **Message gateway reasoning** - Trace how messages, tasks, and agent state should flow between the harness, tools, subagents, and user-facing outputs.
7. **Validation and troubleshooting** - Prefer minimal reproducible runs, logs, configuration checks, and repository tests before recommending architectural changes.

## Inputs / Outputs

**Inputs**

- User objective, target workflow, or Deer Flow issue to diagnose.
- Local repository path or GitHub URL for `bytedance/deer-flow`.
- Available model providers, API keys, runtime environment, sandbox constraints, and tool permissions.
- Desired deliverables such as research reports, code changes, generated artifacts, logs, or deployment notes.

**Outputs**

- A concise setup or operations plan grounded in the Deer Flow repository.
- Concrete commands, configuration notes, or file-level change guidance when implementation is requested.
- A workflow map showing agent stages, tools, skills, subagents, sandbox use, and message handoffs.
- Validation evidence from repository commands, tests, logs, or documented behavior.
- Clear risk notes for unsupported assumptions, missing credentials, unavailable sandboxes, or incomplete documentation.

## References

- Source finding: https://github.com/bytedance/deer-flow
