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
name: pydanticai-2026
description: Build schema-first Python agents with typed tools, structured outputs, provider portability, and observability. Use when application code should stay explicit and validated rather than hidden behind framework magic.
measurable_outcome: Specify a typed agent interface, validated tool schema, and observability path for a new or existing Python agent workflow within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# PydanticAI (2026)

## Workflow

1. Read `references/sources.md` and confirm model-provider support, tool patterns, and current UI/event-stream guidance.
2. Define the agent contract first: input model, output model, tool schemas, and validation boundaries.
3. Keep business logic in normal Python functions and expose only narrow agent entry points.
4. Decide whether UI streaming, multi-agent patterns, or Logfire tracing are required from the start.
5. Add a smoke test that checks schema conformance and one provider failover path when portability matters.

## When to Use

- Structured outputs and type guarantees are central.
- Python service code must remain readable and testable.
- Provider portability matters.
- You want agent behavior that integrates naturally with regular application code and validation.

## Output Requirements

- Return the chosen agent input/output schemas.
- State which provider path is primary and whether fallback is required.
- Include one validation or contract test.
- Include one observability path, preferably Logfire when available.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
