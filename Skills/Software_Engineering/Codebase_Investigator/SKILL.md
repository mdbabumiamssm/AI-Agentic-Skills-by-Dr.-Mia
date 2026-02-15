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
name: codebase-investigator
description: Expertly analyze large codebases to identify patterns, dependencies, and architectural flaws.
keywords:
  - refactoring
  - analysis
  - architecture
  - discovery
  - search
measurable_outcome: Map key components and data flows of a 50k+ LOC repo within 5 minutes, identifying 3+ potential improvements.
license: MIT
metadata:
  author: AI Agentic Skills Team
  version: "2.0.0"
compatibility:
  - system: linux, macos
allowed-tools:
  - list_directory
  - read_file
  - search_file_content
  - glob
  - run_shell_command
---

# Codebase Investigator

A highly specialized skill for understanding complex, unfamiliar software projects. Use this skill to answer "Where is X implemented?" or "How does module A talk to module B?".

## When to Use
- **Onboarding:** Quickly understanding a new project structure.
- **Bug Triaging:** Identifying the root cause location of a bug.
- **Refactoring:** Planning large-scale architectural changes.
- **Documentation:** Generating architectural diagrams or overviews.

## Core Capabilities
1.  **Semantic Search:** Uses `ripgrep` (via `search_file_content`) to find definitions and usages.
2.  **Structure Mapping:** Uses `tree` or `list_directory` to visualize file hierarchy.
3.  **Dependency Analysis:** Inspects `package.json`, `requirements.txt`, etc., to map external libraries.
4.  **Flow Tracing:** Follows function calls across files to understand execution paths.

## Workflow
1.  **Exploration:** Start with `list_directory` and read key files (`README`, `main.py`).
2.  **Targeted Search:** Search for specific keywords related to the query (e.g., "auth", "payment").
3.  **Deep Dive:** Read implementation files of relevant components.
4.  **Synthesis:** Summarize findings into a report or answer.

## Example Usage
```bash
# Agent prompt:
"Investigate how the 'User' model is persisted in the database."
# This triggers a sequence of search_file_content("class User") calls.
```

## Guardrails
- **Read-Only:** Do not modify code during investigation unless explicitly asked.
- **Efficiency:** Use specific search patterns to avoid overwhelming output.
- **Context:** Always check file paths (e.g., `tests/` vs `src/`) to understand context.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->