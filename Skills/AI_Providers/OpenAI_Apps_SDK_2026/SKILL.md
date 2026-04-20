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
name: openai-apps-sdk-2026
description: Build ChatGPT-native apps with the OpenAI Apps SDK, MCP-backed tool servers, widgets, app metadata, and submission-ready integration patterns. Use when the main task is building an app surface inside ChatGPT rather than a raw API-only workflow.
measurable_outcome: Design or update an Apps SDK integration with a clear MCP server plan, widget strategy, metadata posture, and test path within 25 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# OpenAI Apps SDK (2026)

## Workflow

1. Read `references/sources.md` and confirm the current Apps SDK lifecycle before implementation.
2. Define the app around MCP primitives first: tools, schemas, resources, widget metadata, and host interactions.
3. Separate planning from building: use cases, tool definitions, UI components, auth, state, and deployment should be explicit before coding.
4. Use the official examples repository as the starting point for widget and MCP-server composition rather than inventing private conventions.
5. Validate with one local ChatGPT connection or equivalent integration test path before broader rollout.

## When to Use

- You need a ChatGPT-native app rather than only an API workflow.
- The integration requires widgets, tool results with metadata, or conversational UI inside ChatGPT.
- You need an MCP-backed app with submission-ready metadata and security review considerations.
- You need to decide whether a feature belongs in Apps SDK, Codex plugins, or a plain API app.

## Output Requirements

- Return the selected app architecture with a one-line rationale.
- State the MCP server, widget, and metadata strategy.
- State the auth, privacy, and submission assumptions.
- Include a minimal smoke-test plan covering tool listing, tool call, and widget rendering.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
