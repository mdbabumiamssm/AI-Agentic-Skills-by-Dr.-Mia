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
name: 'claude-code-skill-marketplace-ops'
description: 'Assess and safely adopt Claude Code skills/plugins from large marketplaces such as ccpi/tonsofskills, covering discovery, trust, dependencies, conflicts, and installation hygiene.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Claude Code Skill Marketplace Operations

## Overview

This skill guides safe discovery, evaluation, and adoption of third-party Claude Code plugins, skills, and agents from large marketplace repositories. It is intended for agent teams that need operational hygiene before importing broad catalogs into production or shared development environments.

Use it to turn marketplace browsing into a repeatable review workflow: identify candidates, score trust signals, inspect dependencies and instructions, detect namespace conflicts, install in isolation, and document adoption decisions.

## When to Use This Skill

- Evaluating a Claude Code plugin, skill, or agent before installation.
- Reviewing large marketplace catalogs such as ccpi/tonsofskills for safe adoption.
- Comparing multiple third-party skills that overlap with existing internal capabilities.
- Checking dependency, tool, MCP, or shell-command risk in contributed skill packages.
- Investigating namespace collisions, duplicate skill names, or conflicting trigger descriptions.
- Preparing a controlled rollout plan for shared agent workstations or team skill registries.
- Auditing already installed third-party skills for provenance, trust, or maintenance concerns.

## Core Capabilities

1. **Marketplace discovery**
   Search marketplace metadata, repository indexes, package descriptions, topics, and README files to identify candidate Claude Code skills, plugins, and agents that match the requested task.

2. **ccpi/tonsofskills large-catalog operations**
   Review ccpi-managed tonsofskills marketplace catalogs, including the reported 425 plugins, 2,810 skills, and 200 agents, using discovery heuristics, trust scoring, dependency and prompt-injection review, conflict detection, install hygiene, and rollback planning before adoption.

3. **ctx-scale catalog review**
   Evaluate large skill, MCP, and harness catalogs such as ctx with capped recommendations, provenance checks, duplicate detection, dependency and execution-surface review, and conflict scoring before importing catalog-derived skills into curated repositories.

4. **Official Claude plugin directory adoption**
   Treat Anthropic's official claude-plugins-official directory as the preferred trust tier for Claude Code plugin discovery when it satisfies the use case. When an official plugin overlaps with a community marketplace option, compare provenance, capability scope, plugin manifest, dependency footprint, requested permissions, MCP surfaces, and maintenance signals before selecting the install source. Inspect manifests and dependencies before adoption, test selected plugins in isolated repos when possible, pin the selected version or commit before installation when supported, track update hygiene, and document rollback paths for every installed plugin.

- **Official directory trust decisions**
  Prefer Anthropic-managed provenance when available, pin plugin versions or commits, inspect permissions and dependencies, test in a sandbox, document the trust decision, and compare official plugins against community marketplace alternatives before rollout.

5. **Trust and provenance review**
   Record source URL, maintainer identity, repository age, release cadence, stars, issues, license, commit activity, and whether installation instructions are traceable to the upstream project.

6. **Dependency and tool-risk inspection**
   Review package files, scripts, MCP declarations, shell commands, network calls, and requested tool permissions before installation. Flag broad filesystem access, credential handling, destructive commands, hidden downloads, or unpinned dependencies.

7. **Namespace and trigger conflict analysis**
   Compare candidate names, descriptions, commands, MCP server IDs, and skill trigger wording against the local skill/plugin registry to prevent ambiguous activation or accidental override of trusted internal workflows.

8. **Installation hygiene**
   Prefer isolated evaluation first. Install into a scratch or staging environment when possible, avoid bulk importing entire catalogs, pin versions or commit SHAs when supported, and keep a rollback path for every adopted package.

9. **Adoption decision record**
   Produce a concise recommendation: adopt, adopt with restrictions, defer, or reject. Include rationale, risks, required mitigations, owner-reviewed files, and any post-install validation steps.

## Inputs / Outputs

**Inputs**

- Marketplace package name, repository URL, catalog entry, or candidate list.
- Local skill/plugin inventory, if available.
- Installation target: personal sandbox, team workstation, CI image, or shared agent registry.
- Risk tolerance and any organization-specific allowlist, denylist, or required review policy.

**Outputs**

- Candidate summary with source, purpose, scope, and installation mechanism.
- Trust review covering provenance, maintenance signals, license, and repository health.
- Dependency and permission notes, including risky scripts or requested tools.
- Namespace conflict report against known local skills, plugins, commands, and MCP servers.
- Safe installation plan with staging steps, pinning guidance, and rollback notes.
- Final adoption recommendation with required mitigations and validation checks.

## Recommended Workflow

1. **Collect candidate metadata**
   Fetch or read the marketplace entry, upstream repository, package manifest, README, installation instructions, and any included skill/plugin definitions.

2. **Inspect before installing**
   Read package files directly. Do not rely only on marketplace summaries. Pay particular attention to shell scripts, post-install hooks, MCP server configuration, environment variables, and tool permission requests.

3. **Score operational risk**
   Classify risk as low, medium, or high based on provenance, dependency behavior, privilege requests, install scope, namespace overlap, and whether the package executes code during setup.

4. **Check conflicts**
   Compare names, descriptions, commands, MCP server identifiers, and trigger phrases with the local registry. Rename, namespace, or reject candidates that create ambiguous activation.

5. **Stage and validate**
   Install only the selected package into an isolated location when possible. Run a minimal task that exercises the advertised behavior, then inspect changed files and generated configuration.

6. **Record the decision**
   Save the recommendation and review notes in the caller's requested format. Include enough detail for another operator to repeat or reverse the decision.

## Review Checklist

- Source URL and package identity are recorded.
- License and upstream ownership are clear enough for the intended use.
- Installation path, changed files, and rollback action are known.
- Dependencies are pinned or otherwise acceptable for the environment.
- No unexpected credential access, telemetry, hidden downloads, or destructive commands are present.
- Requested tools and MCP servers are necessary for the advertised capability.
- Skill names, plugin IDs, commands, and trigger descriptions do not conflict with trusted local assets.
- A small validation task succeeds before team-wide adoption.
- The final recommendation states adopt, restrict, defer, or reject.

## References

- Source finding: https://github.com/jeremylongshore/claude-code-plugins-plus-skills
- Source finding: https://github.com/anthropics/claude-plugins-official
- Source finding: https://github.com/stevesolun/ctx
