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
name: 'gitagent-protocol'
description: 'Use GitAgent Protocol to define, inspect, and package framework-agnostic, git-native AI agent specifications for repository-based agent workflows.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# GitAgent Protocol

## Overview

This skill helps create, inspect, and validate repository-native AI agent definitions using the GitAgent Protocol. Use it when an agent needs to be represented as portable project files instead of being tied to a single runtime, IDE, or agent framework.

The workflow emphasizes git-friendly structure, clear agent metadata, reusable instructions, and interoperability across coding-agent harnesses such as Codex, Claude Code, OpenClaw, and related tools.

## When to Use This Skill

- Defining an AI agent that should live directly inside a git repository.
- Converting existing agent instructions, prompts, or skills into a portable GitAgent-style layout.
- Reviewing whether a repository contains enough structured context for an agent runtime to load and execute it.
- Creating framework-agnostic agent packages for use across multiple coding assistants or automation systems.
- Auditing agent definitions for missing metadata, unclear instructions, unsafe tool assumptions, or non-portable dependencies.
- Designing repository conventions for reusable agent skills, commands, tools, and documentation.

## Core Capabilities

1. **Repository-native agent modeling**: Identify the files, metadata, and instruction boundaries needed to describe an agent inside a git repository.
2. **Protocol-oriented packaging**: Organize agent definitions so they can be versioned, reviewed, copied, and evolved through normal git workflows.
3. **Framework-agnostic review**: Check agent instructions for assumptions that bind them too tightly to one assistant, CLI, editor, or orchestration framework.
4. **Skill and tool inventory**: Summarize the agent's expected skills, allowed tools, runtime dependencies, and repository access requirements.
5. **Interoperability mapping**: Translate agent intent into portable concepts that can be adapted for Codex, Claude Code, OpenClaw, or other compatible harnesses.
6. **Validation guidance**: Produce a concise checklist for confirming that an agent definition is readable, complete, deterministic enough to execute, and safe to share.
7. **Git-native specification migration**: Handle framework-agnostic GitAgent Protocol definitions by checking repository packaging conventions, validation boundaries, interoperability limits, and migration needs when converting from tool-specific skill formats.
8. **Framework-agnostic GitAgent specification guidance**: Review repository layout, manifest validation, skill and tool declarations, compatibility checks, versioning, packaging, and migration notes for existing agent repositories.
9. **Protocol-first agent specification pattern**: Prefer GitAgent Protocol definitions when the agent needs portable repository layout conventions, packaging, and portability checks across runtimes; keep runtime-specific configs for harness-only execution details.
10. **Current git-native agent specification workflow**: Define a repository layout for portable agent definitions, package reusable skills as git-versioned project assets, run interoperability checks across target harnesses, and migrate from framework-specific agent configs when portability, reviewability, or shared repository execution is required.
11. **Protocol adoption readiness**: Guide migration from ad hoc agent specs into framework-agnostic, git-native definitions by validating manifests, checking repository packaging, confirming interoperability across target harnesses, and reviewing tool, permission, and dependency declarations for security risks.
12. **OpenGAP GitAgent Protocol definition audit**: Check current git-native agent repos for the required `agent.yaml` manifest and `SOUL.md` identity files, optional `RULES.md`, `DUTIES.md`, `AGENTS.md`, `skills/`, `tools/`, `workflows/`, `knowledge/`, `memory/`, `hooks/`, `config/`, `compliance/`, `agents/`, `examples/`, and gitignored `.gitagent/`; verify manifest fields such as `name`, `version`, `description`, `spec_version`, `skills`, `tools`, `extends`, and `dependencies`; bind skills by matching `agent.yaml` skill entries to `skills/<skill-name>/SKILL.md` and tool entries to `tools/<name>.yaml`; run `opengap validate` plus export/import adapter checks for target runtimes; package through normal git review, tagging, and optional `opengap install` dependency resolution; migrate framework-specific configs by moving portable identity, prompts, rules, tool schemas, role/SOD policy, and model preferences into GitAgent files while leaving runtime orchestration, live tool execution, memory I/O, and iterative loops in the source framework.
13. **GitAgent Protocol standardization patterns**: Apply the open-gitagent `gitagent-protocol` finding as a framework-agnostic, git-native standard for documenting repository layout, manifest validation, skill and tool packaging, interoperability checks, and migration notes for teams standardizing agents across clients.
14. **MCP-compatible GitAgent package review**: Use the `open-gitagent/gitagent-protocol` repository as the current framework-agnostic, git-native source when checking repository layout, manifest validation, declared skills and tools, portability across agent harnesses, versioning, and compatibility notes for MCP-backed agent workflows.
15. **Trusted GitAgent Protocol execution review**: Before executing agent definitions from a repository, inspect the git-native layout, validate manifests and declared skills/tools, package only reviewed repository assets, run interoperability checks against the intended agent harnesses, and complete a trust review of permissions, dependencies, setup commands, and runtime assumptions.
16. **open-gitagent specification refresh**: Treat `open-gitagent/gitagent-protocol` as the current framework-agnostic, git-native specification model for defining AI agents; verify repository layout, interoperability checks, git-based packaging guidance, and migration notes from tool-specific skill formats into portable GitAgent files.
17. **Git-native packaging checks**: For existing skill or agent bundles, verify repository layout, manifest validity, framework-agnostic agent definition boundaries, interoperability across target harnesses, and migration notes needed to package them as portable GitAgent Protocol assets.
18. **Discoverable GitAgent package migration**: When migrating client-specific skill formats into GitAgent Protocol, keep the package framework-agnostic and git-native by making agent definitions discoverable in the repository, checking compatibility with target clients, reviewing trust boundaries before execution, and documenting which runtime-specific behavior stays outside the portable package.
19. **Portable GitAgent protocol packaging**: Package agent repositories with a clear git-native layout, validated agent manifest, framework-agnostic skill and tool definitions, interoperability checks for target runtimes, trust review of permissions and dependencies, and migration notes for teams converting existing agent repositories into portable definitions.

## Inputs / Outputs

**Inputs**

- A repository path, GitHub URL, or file tree containing agent definitions, prompts, skills, or automation instructions.
- Optional target runtime context, such as Codex, Claude Code, OpenClaw, or another agent harness.
- Optional requirements for allowed tools, repository permissions, setup commands, environment variables, or validation steps.
- Existing documentation describing the agent's purpose, operating constraints, expected tasks, and success criteria.

**Outputs**

- A GitAgent-oriented agent definition plan or review.
- A recommended repository layout for portable agent files and supporting documentation.
- A gap analysis covering missing metadata, unclear instructions, runtime assumptions, and portability risks.
- Concrete edits or markdown snippets when the user asks to author or revise agent specification files.
- A validation checklist that can be executed manually or through repository automation.

## References

- GitAgent Protocol repository: https://github.com/open-gitagent/gitagent-protocol
- Source finding: https://github.com/open-gitagent/gitagent-protocol
