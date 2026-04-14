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
name: a2a-protocol-2026
description: Design interoperable agent services using the Agent2Agent protocol. Use when independently hosted agents need to discover each other, exchange tasks, and collaborate without sharing a framework runtime.
measurable_outcome: Define an A2A-compatible agent surface, task exchange pattern, and security posture for a concrete workflow within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# A2A Protocol (2026)

## Workflow

1. Read `references/sources.md` and confirm the current A2A concepts, enterprise requirements, and official SDK paths.
2. Define the agent card and public capabilities before discussing implementation details.
3. Decide the task lifecycle, attachment strategy, and status updates that remote agents must exchange.
4. Keep security and network boundaries explicit: HTTPS, auth, and opaque remote-agent assumptions.
5. Add an interop smoke test between two simple agents before layering provider-specific logic on top.

## When to Use

- Remote agents must communicate across framework boundaries.
- A client should be able to discover and invoke another agent as a service.
- You want a standard protocol instead of a custom RPC contract.
- The design must keep internal memory and tool state private while still enabling collaboration.

## Output Requirements

- Return the agent card fields that must be exposed.
- State the task exchange pattern and status model.
- Include one transport or security assumption.
- Include one note on how A2A will coexist with MCP, if both are present.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
