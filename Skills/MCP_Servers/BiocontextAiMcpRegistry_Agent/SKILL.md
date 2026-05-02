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
name: 'biocontext-ai-mcp-registry'
description: 'Discover, evaluate, and integrate biomedical MCP servers from the BioContextAI Registry — a curated catalog of Model Context Protocol servers for bioinformatics, systems biology, and biomedical agentic AI.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# BioContextAI MCP Registry

## Overview

This skill helps an agent browse, vet, and wire up biomedical Model Context Protocol (MCP) servers listed in the **BioContextAI Registry** — a domain-specific catalog focused on bioinformatics, systems biology, and biomedical AI tooling. It complements general MCP operations skills by adding a curated discovery layer so agents can stand up biomedical tool stacks (e.g., variant lookup, pathway enrichment, literature retrieval) without hand-rolling each integration.

## When to Use This Skill

- Building a biomedical or life-sciences agent that needs structured access to bio-databases, ontologies, or analysis tools through MCP.
- Selecting an MCP server for a specific biomedical task (sequence analysis, drug-target lookup, single-cell tools, clinical evidence retrieval).
- Auditing an existing biomedical agent stack to replace bespoke API wrappers with community-maintained MCP servers.
- Contributing a new MCP server back to the registry and needing to follow its submission conventions.
- Cross-checking that a candidate MCP server is actually maintained, scoped to biomedical use, and safe to call from an agent.

## Core Capabilities

1. **Registry discovery** — Fetch and parse the BioContextAI Registry index from GitHub to enumerate available biomedical MCP servers, their domains, and maintainers.
2. **Server triage** — For a given biomedical task, shortlist candidate MCP servers from the registry by topic (e.g., genomics, proteomics, pharmacology, literature), language, and activity signals (stars, recent commits).
3. **Capability extraction** — Read each candidate server's README/manifest to enumerate exposed MCP tools, resources, and prompts, and summarize their input/output contracts.
4. **Stack assembly** — Produce a minimal client configuration (JSON snippet for Claude Desktop / Claude Code / generic MCP host) that wires the selected servers into the agent runtime.
5. **Compliance & safety check** — Flag servers that proxy clinical/PHI data, call paid APIs, or require credentials so the user can review licensing, rate-limit, and privacy implications before enabling them.
6. **Submission scaffolding** — When the user has built a new biomedical MCP server, generate the metadata entry needed to submit it to the registry (name, description, tags, repo URL, install command).

## Inputs / Outputs

**Inputs**
- A biomedical task description or research question (e.g., "look up gene–disease associations and recent PubMed evidence").
- Optional constraints: preferred language runtime (Python/Node/Rust), local-only vs. hosted servers, license requirements, credential availability.
- Optional path to an existing MCP host configuration file to extend.

**Outputs**
- A ranked shortlist of registry servers with one-line capability summaries and source repo links.
- An MCP host config snippet (JSON) ready to drop into the user's client (Claude Desktop `claude_desktop_config.json`, Claude Code `.mcp.json`, or equivalent).
- A compliance note listing credential, rate-limit, and data-sensitivity considerations per selected server.
- (Optional) A registry submission stub when contributing a new server.

## References

- Source finding: BioContextAI Registry — https://github.com/biocontext-ai/registry
- Model Context Protocol specification — https://modelcontextprotocol.io
- MCP reference servers — https://github.com/modelcontextprotocol/servers
- Anthropic MCP announcement — https://www.anthropic.com/news/model-context-protocol
