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
7. Vetted-directory review: Use the official GitHub directory metadata, plugin trust signals, install hygiene checks, compatibility testing, and migration guidance to move from unofficial marketplaces or community catalogs toward official Anthropic-managed Claude Code plugins.
8. Official directory adoption: Prioritize vetted plugin discovery from `anthropics/claude-plugins-official`; record provenance checks, dependency review, installation hygiene steps, version tracking, and migration notes when replacing community plugin catalogs.
9. Preferred vetted source workflow: Treat Anthropic's official `claude-plugins-official` GitHub directory as the preferred plugin discovery source, then document trust review, dependency hygiene, installation boundaries, and migration guidance before using unofficial marketplaces.
10. Catalog-governed operations: Vet candidate plugins against the managed catalog, compare any community marketplace alternatives, review install permissions and bundled MCP servers, pin trusted versions where possible, and document fallback paths for deprecated plugins.
11. Official marketplace migration: Start plugin discovery from the official Anthropic-managed `anthropics/claude-plugins-official` directory, inspect dependencies and plugin manifests before install, pin approved versions or commits, keep installation state clean and reversible, and document migration decisions when replacing community marketplace plugins.
12. Preferred trust baseline: Use the official Anthropic-managed `anthropics/claude-plugins-official` directory as the baseline for plugin discovery, vetting, dependency review, installation hygiene, and migration away from uncurated marketplace entries.
13. Current official-plugin operations: Vet official plugin entries against their repository and manifest before adoption, inspect permissions, dependencies, and bundled MCP/tool declarations, pin approved versions or commits, test installs in isolated workspaces, handle conflicts with existing skills or MCP tools, and plan staged migrations from unofficial plugin sources.
14. Official plugin directory operations: Prefer Anthropic-managed sources when available, verify plugin provenance, inspect permissions and dependencies, pin approved versions or commits, document migration from community plugins, and keep installation approval explicit.
15. Vetted listing lifecycle: Use `anthropics/claude-plugins-official` as the vetted discovery source, capturing the listing URL, visible repository metadata, trust review, dependency checks, install hygiene decisions, approved version or commit, and migration plan from unofficial marketplaces.
16. Directory-vetting adoption checklist: Before adopting a plugin from `anthropics/claude-plugins-official`, verify the trusted source, review plugin scope against the requested workflow, inspect dependency and permission hygiene, check visible update cadence, test for conflicts with existing skills, plugins, or MCP tools, and define rollback steps.
17. Authoritative directory preference: Treat `anthropics/claude-plugins-official` as the authoritative vetted plugin directory when a matching official plugin exists; prefer it over marketplace packages when provenance is clearer, dependencies and install steps can be reviewed, migration impact is documented, and the official option covers the required workflow.
18. Anthropic-managed official plugin workflow: Use the GitHub `anthropics/claude-plugins-official` directory as a vetted source for Claude Code plugin discovery; before installation, perform trust review, inspect dependencies and permission requests, keep installs clean and reversible, avoid conflicts with existing plugins, skills, or MCP servers, and document migration decisions when replacing informal plugin directories.
19. Trusted adoption baseline: Treat the official Anthropic-managed `anthropics/claude-plugins-official` directory as the baseline for Claude Code plugin adoption; distinguish official entries from marketplace plugins, inspect dependencies and permission requests, pin approved versions or commits, avoid namespace conflicts with existing plugins, skills, and MCP servers, and plan migrations from unofficial plugin bundles.

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
- https://github.com/anthropics/claude-plugins-official
- https://github.com/anthropics/claude-plugins-official
- GitHub finding, published 2026-05-04 by `anthropics`: https://github.com/anthropics/claude-plugins-official
- Official Anthropic-managed Claude Code plugin directory: https://github.com/anthropics/claude-plugins-official
