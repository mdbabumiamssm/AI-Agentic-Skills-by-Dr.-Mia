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
name: 'officecli-agent-office-automation'
description: 'Use OfficeCLI to read, edit, and automate Word, Excel, and PowerPoint files through a single agent-friendly CLI without requiring Microsoft Office.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# OfficeCLI Agent Office Automation

## Overview

This skill guides agents through using OfficeCLI as a unified command-line workflow for Word, Excel, and PowerPoint deliverables. It is useful when a task needs office document inspection, editing, or automation without relying on a local Microsoft Office installation.

OfficeCLI is especially relevant for mixed Office workflows because it is described as an agent-native, open-source, single-binary tool for reading, editing, and automating Office files.

## When to Use This Skill

- A user asks to inspect, summarize, modify, or automate `.docx`, `.xlsx`, or `.pptx` files.
- The task spans multiple Office formats and a unified CLI workflow is preferable to separate format-specific tools.
- The environment does not have Microsoft Office installed, or the user wants a reproducible command-line path.
- The user needs batch processing, validation, or scripted transformations across Office deliverables.
- Existing `docx`, `xlsx`, or `pptx` skills are too narrow for a mixed document bundle.

## Core Capabilities

1. **Office file discovery and triage**: Identify target Word, Excel, and PowerPoint files, confirm paths, and determine whether the requested work is read-only, editing, or batch automation.
2. **Runtime command discovery**: Check whether `officecli` is installed, run its help output, and use upstream documentation when needed instead of assuming unsupported subcommands.
3. **Single-binary Office automation**: Use OfficeCLI's no-Office-install dependency model for agent-safe Word, Excel, and PowerPoint read/edit flows, verifying outputs by re-reading or inspecting the generated files before handoff.
4. **Agent-friendly read-edit-write loop**: Locate the `officecli` binary before use, inspect available commands for Word, Excel, and PowerPoint operations, assume no Microsoft Office runtime is required, write edits to copies or an output directory by default, and re-read outputs through OfficeCLI before handoff.
5. **Unified document reading**: Extract or inspect document content, workbook structure, sheets, slides, metadata, and other CLI-supported representations for agent analysis.
6. **CLI-supported editing**: Apply OfficeCLI-supported edits to documents while preserving the original inputs unless the user explicitly authorizes in-place mutation.
7. **Mixed-format automation**: Coordinate Word, Excel, and PowerPoint operations in one repeatable shell workflow for reports, analyses, decks, and accompanying spreadsheets.
8. **Safe batch editing**: For multi-file operations, enumerate targets first, write to copies or a dedicated output directory by default, avoid in-place mutation without explicit authorization, and keep per-file command results for rollback or review.
9. **Fidelity validation and handoff**: Re-read or inspect outputs after edits, confirm files remain openable and structurally intact, compare key document features such as Word text coverage, Excel sheets/ranges, or PowerPoint slide counts against the intended changes, capture command failures, and return produced file paths plus a concise summary of changes.

## Inputs / Outputs

**Inputs**

- User request describing the desired Office document operation.
- One or more `.docx`, `.xlsx`, or `.pptx` paths.
- Optional constraints such as target filenames, output directory, formatting requirements, sheet or slide names, and whether edits may be in place.
- Optional source data, replacement text, tables, charts, or instructions to merge into the Office files.

**Outputs**

- Parsed document content, structure summaries, extracted text, or metadata when the task is read-only.
- Edited Word, Excel, or PowerPoint files when the task requests modification.
- Batch automation results such as generated deliverables, transformed copies, or command logs when useful for verification.
- A final response listing changed files, important assumptions, and any OfficeCLI limitations or command failures encountered.

## References

- OfficeCLI GitHub repository: https://github.com/iOfficeAI/OfficeCLI
