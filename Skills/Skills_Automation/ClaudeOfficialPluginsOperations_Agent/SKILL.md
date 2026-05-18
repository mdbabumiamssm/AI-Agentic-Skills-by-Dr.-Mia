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
name: 'claude-official-plugins-operations'
description: 'Use Anthropic-managed claude-plugins-official for vetted Claude Code plugin discovery, trust review, installation hygiene, dependency assessment, and migration planning.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Claude Official Plugins Operations

## Overview
This skill guides work that uses the official Anthropic-managed `claude-plugins-official` directory as the primary source for Claude Code plugin discovery and operational review. It helps distinguish official-source guidance from community catalogs while supporting safe plugin selection, installation hygiene, dependency assessment, and migration planning in agentic skill/plugin ecosystems.

## When to Use This Skill
- A user asks which Claude Code plugins to install, compare, trust, audit, update, or migrate.
- A workflow needs official Anthropic-managed plugin information before considering community plugins.
- A plugin adoption decision requires provenance, dependency, permission, or maintenance review.
- A team needs a safe install, upgrade, rollback, or validation plan for Claude Code plugins.
- An existing skill/plugin ecosystem needs migration planning toward official Claude Code plugins.
- A user asks to separate official directory guidance from third-party marketplace or repository claims.

## Core Capabilities
1. Official-source discovery: Query or inspect the `anthropics/claude-plugins-official` directory first, recording plugin names, repository links, stated purpose, and source context.
2. Trust and provenance review: Confirm whether plugin information comes from the official Anthropic-managed directory, then flag any claims that come from external or community sources.
3. Dependency assessment: Review plugin manifests, package files, scripts, MCP server declarations, and installation instructions for runtime dependencies, network access, credential needs, and local system requirements.
4. Installation hygiene: Produce conservative preflight steps, backup guidance, install commands, environment-variable handling, validation checks, and rollback notes without assuming plugins are safe merely because they are discoverable.
5. Migration planning: Map current skills, MCP servers, or plugins to official Claude Code plugin options, identifying overlaps, gaps, compatibility risks, and staged adoption steps.
6. Operational reporting: Return concise findings that separate evidence, recommendation, residual risk, and next actions so users can make adoption decisions quickly.

## Inputs / Outputs
Inputs:
- User goal, target workflow, or plugin capability request.
- Candidate plugin names, repository URLs, directory entries, or current local plugin configuration.
- Target Claude Code environment, operating system constraints, package manager constraints, and credential policies when available.
- Existing skills, MCP servers, automation scripts, or plugin inventory when migration is requested.

Outputs:
- Official-source plugin shortlist with links and the reason each candidate matches the user goal.
- Trust review covering provenance, dependency surface, permission needs, maintenance signals visible from the source, and notable unknowns.
- Installation or update plan with preflight checks, install steps, validation commands, backup points, and rollback guidance.
- Migration plan mapping current capabilities to official plugin options, including retained gaps and staged rollout steps.
- Clear distinction between official directory evidence and any supplemental third-party information used for context.

## References
- Source finding: `anthropics/claude-plugins-official`, official Anthropic-managed directory of high quality Claude Code Plugins: https://github.com/anthropics/claude-plugins-official
