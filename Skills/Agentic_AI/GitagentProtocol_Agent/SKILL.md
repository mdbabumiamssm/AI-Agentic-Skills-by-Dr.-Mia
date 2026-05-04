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
