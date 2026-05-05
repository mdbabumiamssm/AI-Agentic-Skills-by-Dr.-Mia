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
name: 'open-design-agentic-prototyping'
description: 'Use Open Design to plan safe local-first agentic workflows for design systems, prototypes, previews, and HTML/PDF/PPTX/MP4 exports.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Open Design Agentic Prototyping

## Overview

This skill guides safe evaluation and workflow integration of Open Design, a local-first open-source design and prototyping runtime for agentic coding environments. Use it to translate product or design intent into inspectable web, desktop, mobile, slide, image, video, or HyperFrame prototype workflows while preserving local control, sandboxed preview, and explicit export expectations.

## When to Use This Skill

- The user wants a local-first alternative to hosted AI design tools.
- The user needs an agent-assisted prototype for web, desktop, mobile, slides, images, videos, or HyperFrames.
- The task involves evaluating or adopting Open Design inside Codex, Claude Code, Cursor, Gemini CLI, OpenCode, Qwen, Copilot, Hermes, or Kimi CLI workflows.
- The user asks for design-system-driven generation, brand-grade UI prototyping, or reusable design workflows.
- The user needs sandboxed preview or export planning for HTML, PDF, PPTX, or MP4 deliverables.
- The user wants a safe adoption checklist before running a broad third-party agentic design tool locally.

## Core Capabilities

1. **Repository orientation**: Inspect the Open Design repository, current setup guidance, supported runtimes, and documented skill surfaces before recommending commands.
2. **Workflow selection**: Map the user's deliverable to an appropriate Open Design path, such as prototype generation, design-system application, preview, or export.
3. **Local-first setup planning**: Identify prerequisites, environment variables, package manager expectations, and BYOK considerations without assuming hosted services.
4. **Sandboxed preview discipline**: Prefer preview flows that isolate generated code and make artifacts inspectable before reuse in a production repository.
5. **Export planning**: Clarify target output format, acceptance criteria, and handoff artifacts for HTML, PDF, PPTX, MP4, or other supported outputs.
6. **Integration safety review**: Check licensing, dependency footprint, generated-code boundaries, secret handling, and repository cleanliness before adoption.

## Inputs / Outputs

**Inputs**

- User goal, target platform, desired artifact type, visual direction, and brand or design-system constraints.
- Existing repository context, if the prototype must integrate with an app or design system.
- Runtime preference, such as Codex, Claude Code, Cursor, Gemini CLI, OpenCode, Qwen, Copilot, Hermes, or Kimi CLI.
- Export requirements, including format, dimensions, interaction depth, and review deadline.

**Outputs**

- A concise Open Design workflow plan with setup, generation, preview, validation, and export steps.
- Prototype or design-system integration instructions grounded in the repository documentation.
- A safety checklist covering dependencies, secrets, sandboxing, generated assets, and licensing review.
- Valid deliverable artifacts or clear next commands for producing HTML, PDF, PPTX, MP4, or other documented outputs.

## References

- Source finding: https://github.com/nexu-io/open-design
