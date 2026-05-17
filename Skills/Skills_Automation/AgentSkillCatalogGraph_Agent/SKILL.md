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
name: 'agent-skill-catalog-graph'
description: 'Use graph-backed catalogs such as ctx to discover, compare, and recommend skills, agents, MCPs, and harnesses for LLM automation workflows.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Agent Skill Catalog Graph

## Overview

This skill guides graph-backed discovery and selection of agent skills, MCP servers, agents, and harnesses using catalog projects such as `stevesolun/ctx`. It is useful when a workflow needs more than keyword search: dependency checks, overlap analysis, conflict detection, and recommendation of a small execution-ready toolset from a large catalog.

## When to Use This Skill

- The user asks which skills, MCPs, agents, or harnesses should be used for a task.
- A repository or workspace contains many candidate skills and needs catalog-level triage.
- A task requires detecting duplicated, overlapping, missing, or conflicting skill capabilities.
- The user wants a capped recommendation set rather than a broad marketplace dump.
- An automation workflow needs a repeatable method for selecting LLM tools from graph metadata.
- The user references `ctx`, LLM-wiki graphs, skill catalogs, MCP catalogs, or harness recommendations.

## Core Capabilities

1. **Catalog ingestion** - Locate the relevant catalog source, fetch or read its metadata, and identify available skills, agents, MCPs, harnesses, and graph files.
2. **Requirement modeling** - Convert the user's task into explicit capability requirements, constraints, preferred environments, and exclusion criteria.
3. **Graph-backed discovery** - Use relationships among catalog nodes to find adjacent tools, complementary skills, dependency paths, and likely substitutions.
4. **Dependency and conflict analysis** - Check whether recommended items require missing runtimes, incompatible tools, duplicate responsibilities, or mutually exclusive assumptions.
5. **Recommendation capping** - Return a concise set of candidates ranked by task fit, coverage, operational risk, and setup burden instead of listing every match.
6. **Workflow handoff** - Produce practical next steps: install/use candidates, verify prerequisites, run sample commands, or document why an item was rejected.

## Inputs / Outputs

**Inputs**

- User goal, task description, or target automation workflow.
- Optional existing skill, MCP, agent, or harness inventory from local files or a remote catalog.
- Optional constraints such as allowed tools, runtime, security posture, offline use, language, repository style, or maximum number of recommendations.
- Optional source reference, such as the `stevesolun/ctx` GitHub repository or a local clone of that catalog.

**Outputs**

- A ranked shortlist of recommended skills, MCPs, agents, or harnesses.
- A short rationale for each recommendation, including the capability covered.
- Any detected dependencies, conflicts, duplicates, or missing prerequisites.
- A capped execution plan that states what to use first, what to defer, and what to avoid.
- Optional shell commands or file paths for inspecting a local catalog when the user permits command execution.

## References

- Source finding: https://github.com/stevesolun/ctx
