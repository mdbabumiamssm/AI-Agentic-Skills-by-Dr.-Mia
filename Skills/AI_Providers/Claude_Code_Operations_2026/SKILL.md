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
name: claude-code-operations-2026
description: Operate Anthropic Claude Code across terminal use, agent teams, Agent SDK, GitHub Actions, plugins, and CI-based repository automation. Use when the main task is coding-agent execution, developer workflow acceleration, or repository-native automation with Anthropic tooling.
measurable_outcome: Select and configure a Claude Code operating path with explicit install, instruction-file, and automation choices within 20 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Claude Code Operations (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Claude Code install and automation surfaces before coding.
2. Choose whether the task belongs in the interactive terminal, agent teams, the Agent SDK, or GitHub Actions.
3. Define repository-local behavior control with `CLAUDE.md`, tool permissions, and CI permissions before enabling broad code changes.
4. Prefer the official GitHub Action for repository automation and the Agent SDK only when you need custom orchestration beyond the built-in product flow.
5. Validate with one local coding task or one CI-triggered repository task.

## Core Capabilities

- Use Anthropic's official Claude Code plugins directory as the preferred Anthropic-managed source for vetted plugin discovery; verify plugin manifests, requested permissions, MCP server exposure, skills interaction, and automation side effects before installing; install only from trusted official sources, pin plugin versions where reproducibility matters, keep official plugins separate from community marketplace assets in organization allowlists and documentation, and migrate away from unverified marketplaces when an official plugin covers the need.

## When to Use

- You need Anthropic's coding agent rather than the general Claude API.
- The workflow is terminal-native and repository-aware.
- You need GitHub-native automation, issue triage, or PR implementation using official Claude Code surfaces.
- You need Anthropic-managed agent teams or Agent SDK control for coding tasks.

## Output Requirements

- Return the chosen Claude Code surface with a one-line rationale.
- State the required install or execution path.
- State the repo-local instruction and permission boundary.
- Include a minimal smoke-test plan covering one interactive or CI-triggered task.

## References

- https://github.com/anthropics/claude-plugins-official


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
