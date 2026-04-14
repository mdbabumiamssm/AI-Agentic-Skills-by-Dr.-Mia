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
name: mcp-registry-2026
description: Design and curate Model Context Protocol tool surfaces, server registries, and shared resource access. Use when agents need portable tool/resource interoperability across clients, runtimes, or providers.
measurable_outcome: Define an MCP server strategy, transport choice, and validation path for a concrete tool surface within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# MCP Registry (2026)

## Workflow

1. Read `references/sources.md` and confirm the current MCP concepts, transports, registry surface, and server patterns.
2. Decide whether the need is server discovery, new server authoring, local stdio integration, or remote HTTP deployment.
3. Define the tool and resource contracts before wiring them into any one framework.
4. Record authentication, approval, and network-boundary assumptions explicitly for every server surface.
5. Add a smoke test that checks discovery, one tool invocation, and one failure mode such as auth rejection or transport disconnect.

## When to Use

- Multiple agent clients should share the same tools or resources.
- Tool access must remain portable across providers or frameworks.
- You need a principled boundary between agent orchestration and tool execution.
- A local or remote MCP server needs to be selected, registered, or exposed safely.

## Output Requirements

- Return the MCP transport choice and why.
- State whether the server is local, remote, or mixed.
- Include one auth or approval assumption.
- Include one validation step for registry discovery or tool invocation.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
