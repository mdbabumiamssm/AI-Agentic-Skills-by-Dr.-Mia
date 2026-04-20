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

---
name: openai-codex-operations-2026
description: Use OpenAI Codex as a coding-agent surface across the app, CLI, IDE extension, SDK, MCP, plugins, skills, subagents, and GitHub automation. Use when the main task is repository automation, code understanding, code changes, or developer workflow acceleration rather than a generic API integration.
measurable_outcome: Select and configure a Codex operating surface with clear security assumptions, repository instruction files, and a minimal smoke-test plan within 20 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# OpenAI Codex Operations (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Codex product surface before implementing anything.
2. Choose the execution surface first: app, CLI, IDE extension, SDK, or GitHub Action.
3. Define repository controls before broad autonomy: `AGENTS.md`, approvals, MCP servers, plugins, skills, and subagents.
4. Use Codex for codebase-local work and direct OpenAI API or Agents SDK primitives when you need lower-level custom orchestration.
5. Validate with one minimal task that exercises the intended surface and permission model.

## When to Use

- You need a first-class coding agent rather than a raw text-generation workflow.
- The work involves real repositories, commands, patches, reviews, or PR automation.
- You need repo-local behavior control through `AGENTS.md`, rules, hooks, skills, MCP, plugins, or subagents.
- You need to decide between local Codex usage and cloud-delegated Codex work.

## Output Requirements

- Return the selected Codex surface with a one-line rationale.
- State the required repository control files or runtime configuration.
- State the approval and security boundary for code execution.
- Include a minimal smoke-test plan covering one read task and one write task.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
