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
name: smolagents-2026
description: Build lightweight code-first or tool-calling agents with Hugging Face smolagents. Use when minimal abstractions, code agents, sandboxed execution, and model flexibility matter more than a heavy orchestration runtime.
measurable_outcome: Select a smolagents agent type, sandbox boundary, and model path for a concrete workload within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# smolagents (2026)

## Workflow

1. Read `references/sources.md` and confirm the current smolagents support for code agents, tool-calling agents, MCP, and secure execution.
2. Choose the right agent type first: `CodeAgent` for code-native action loops or `ToolCallingAgent` for standard tool-calling.
3. Define the execution boundary before enabling autonomy: local Python, Docker, or an external sandbox provider.
4. Decide which model path is primary: Hugging Face, LiteLLM-backed providers, or local Transformers/Ollama.
5. Add a smoke test that covers one agent step, one tool or code path, and one sandbox or safety assumption.

## When to Use

- You want a minimal Python framework for multi-step agents.
- Code agents are a better fit than large orchestration graphs.
- Model portability matters across Hugging Face, provider APIs, and local runtimes.
- Secure execution boundaries must be explicit when the agent writes or runs code.

## Output Requirements

- Return the selected smolagents agent type and why.
- State the model path and sandbox boundary.
- Include one memory, telemetry, or debugging note.
- Include one safety guardrail for code execution or tool exposure.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
