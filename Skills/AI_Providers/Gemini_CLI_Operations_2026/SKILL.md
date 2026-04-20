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
name: gemini-cli-operations-2026
description: Operate Google's Gemini CLI for terminal-native coding, MCP integration, grounding, checkpointing, and headless automation. Use when the main task is developer workflow acceleration through Gemini CLI rather than a general Gemini API integration.
measurable_outcome: Select and configure a Gemini CLI execution mode with explicit auth, context-file, and automation assumptions within 20 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Gemini CLI Operations (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Gemini CLI release and documentation state before coding.
2. Choose the execution mode first: interactive terminal, headless scripting, or GitHub workflow integration.
3. Define auth and behavior controls explicitly, including Google sign-in, API-key paths, `GEMINI.md`, MCP servers, and trusted-folder assumptions.
4. Use grounding, checkpointing, or GitHub integration only where they materially improve the workflow rather than by default.
5. Validate with one small repository task in the selected mode.

## When to Use

- You need a Google-backed terminal coding agent rather than the generic Gemini API.
- The workflow benefits from MCP, grounding, checkpointing, or `GEMINI.md` project context.
- You need headless CLI automation or GitHub-native help from Gemini CLI.
- You want a coding-agent surface that stays distinct from ADK-based application runtimes.

## Output Requirements

- Return the chosen Gemini CLI mode with a one-line rationale.
- State the auth path and project-context file assumptions.
- State whether MCP, grounding, or checkpointing is enabled and why.
- Include a minimal smoke-test plan covering one prompt and one repository action.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
