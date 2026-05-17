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
name: 'claude-code-delivery-lifecycle-skills'
description: 'Operate the claude-code-skills delivery lifecycle suite for agile workflows, AI review, bootstrapping, docs, audits, optimization, MCP editing, graphs, and SSH.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Claude Code Delivery Lifecycle Skills

## Overview

Use this skill to apply the `claude-code-skills` delivery-lifecycle stack to software projects that need structured agentic execution from bootstrap through review, documentation, audit, and optimization. It focuses on safely orchestrating bundled Claude Code skills and MCP servers while preserving clear operational boundaries for repository changes, remote access, and verified edits.

This skill is most useful when a project needs a repeatable engineering workflow rather than a single ad hoc code change.

## When to Use This Skill

- The user asks to use or evaluate `levnikolaevich/claude-code-skills`.
- A task requires a full delivery lifecycle: backlog grooming, planning, implementation support, review, documentation, and release-oriented handoff.
- A project needs bootstrap assistance, codebase audit, performance optimization, or documentation generation through Claude Code workflows.
- The work benefits from multi-model AI review or structured agent review gates.
- The repository needs hash-verified line editing, code knowledge graph analysis, or remote SSH workflows through bundled MCP servers.
- The user asks for operational boundaries around Claude Code plugins, skills, MCP servers, or delivery automation.

## Core Capabilities

1. **Lifecycle orchestration** - Select the appropriate delivery phase and run a bounded workflow for planning, implementation support, review, documentation, audit, or optimization.
2. **Agile pipeline support** - Convert user goals into scoped work items, acceptance criteria, delivery steps, and review checkpoints.
3. **Multi-model review coordination** - Use AI review as a quality gate for implementation plans, diffs, architecture concerns, security risks, and regression risks.
4. **Project bootstrap** - Initialize or assess repository structure, development conventions, documentation expectations, and early automation requirements.
5. **Documentation generation** - Produce repository-grounded documentation such as onboarding notes, architecture summaries, operational guides, and change handoffs.
6. **Codebase audit** - Inspect project structure, dependency surfaces, maintainability risks, architecture drift, and likely hotspots before recommending changes.
7. **Performance optimization workflow** - Identify measurable bottlenecks, choose profiling or benchmark checks, make scoped improvements, and record before/after evidence when available.
8. **Verified editing boundary** - Prefer hash-verified or line-verified editing workflows when using `hex-line` style operations so changes apply to the intended source text.
9. **Code graph analysis boundary** - Use `hex-graph` style code knowledge graphs for navigation, dependency discovery, and impact analysis; do not treat graph output as a substitute for reading critical source files.
10. **Remote SSH boundary** - Use `hex-ssh` style remote execution only when explicitly authorized, with narrow commands, clear targets, and no unrelated system changes.

## Inputs / Outputs

**Inputs**

- User goal, repository path, branch or working tree constraints, and delivery phase.
- Existing project files, issue or ticket context, documentation, tests, and build commands.
- Optional access details for plugin skills, MCP servers, AI review tools, code graph tooling, verified editing, or SSH targets.
- Explicit safety constraints such as files to avoid, commands not to run, deployment limits, and review requirements.

**Outputs**

- A scoped workflow plan tied to the requested delivery phase.
- Concrete project artifacts such as implementation notes, review findings, audit reports, documentation drafts, optimization reports, or bootstrap checklists.
- Edited code or documentation only when the user asks for implementation and the repository state permits safe changes.
- Validation evidence from tests, linters, build checks, profiling, code graph inspection, or AI review summaries when available.
- A concise handoff with changed files, unresolved risks, and recommended next steps.

## References

- Source finding: https://github.com/levnikolaevich/claude-code-skills
- Model Context Protocol: https://github.com/modelcontextprotocol
- Tree-sitter parsing toolkit: https://github.com/tree-sitter/tree-sitter
