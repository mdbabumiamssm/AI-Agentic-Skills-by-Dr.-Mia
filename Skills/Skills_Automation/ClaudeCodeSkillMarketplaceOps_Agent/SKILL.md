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

- **ccpi/tonsofskills scale and adoption checks**
  Treat `jeremylongshore/claude-code-plugins-plus-skills`, tonsofskills.com, and the ccpi CLI as a large third-party marketplace source with reported scale signals of 425 plugins, 2,810 skills, 200 agents, and 2,095 GitHub stars. Use those signals for discovery prioritization only; still perform trust review, dependency inspection, namespace conflict handling, version or commit pinning when supported, staged install validation, rollback planning, and update-diff review before adoption.

- **ccpi marketplace adoption gate**
  For ccpi/tonsofskills candidates, score trust before installation, review direct and transitive dependencies, check skill/plugin/command/MCP namespace conflicts, pin versions or commits when supported, install into a staging profile with clear uninstall notes, and fork rather than install directly when local hardening, dependency pinning, namespace edits, or policy-specific restrictions are required.

- **ccpi/tonsofskills collection triage**
  For large third-party ccpi/tonsofskills collections, triage catalog entries before installation, score trust and provenance, review dependencies plus requested tools, permissions, scripts, and MCP surfaces, check skill/plugin/command namespace conflicts, use a lockfile or version/commit pinning workflow when supported, and keep rollback notes that identify install paths, changed files, and uninstall steps.

- **ccpi catalog acceptance and rejection criteria**
  For ccpi/tonsofskills skill packs, combine catalog trust scoring with dependency and permission review, version or commit pinning, rollback notes, and namespace conflict handling before adoption. Reject low-quality or duplicate packs when provenance is unclear, instructions are too broad to audit, dependencies or requested permissions are unjustified, triggers conflict with trusted local skills, or an existing vetted package already covers the same capability.

- **ccpi scale update and removal controls**
  For tonsofskills-scale catalogs installed through ccpi or equivalent package-manager flows, keep a lockfile or reviewed manifest of adopted plugins, skills, agents, versions or commits, dependencies, namespace decisions, install paths, and removal commands. Before updates, diff manifests, scripts, dependencies, and permission surfaces, then rerun trust scoring, dependency review, and conflict detection. Before removal, verify no local skill, plugin, agent, command, or MCP dependency still requires the package, test uninstall in staging, and document shared-registry cleanup.

3. **ctx-scale catalog review**
   Evaluate large skill, MCP, and harness catalogs such as ctx with capped recommendations, provenance checks, duplicate detection, dependency and execution-surface review, and conflict scoring before importing catalog-derived skills into curated repositories.

4. **Official Claude plugin directory adoption**
   Treat Anthropic's official claude-plugins-official directory as the preferred trust tier for Claude Code plugin discovery when it satisfies the use case. When an official plugin overlaps with a community marketplace option, compare provenance, capability scope, plugin manifest, dependency footprint, requested permissions, MCP surfaces, and maintenance signals before selecting the install source. Inspect manifests and dependencies before adoption, test selected plugins in isolated repos when possible, pin the selected version or commit before installation when supported, track update hygiene, and document rollback paths for every installed plugin.

- **Official directory trust decisions**
  Prefer Anthropic-managed provenance when available, pin plugin versions or commits, inspect permissions and dependencies, test in a sandbox, document the trust decision, and compare official plugins against community marketplace alternatives before rollout.

- **Official source-tier separation**
  Treat `anthropics/claude-plugins-official` as a higher-trust discovery tier when a matching plugin exists, while keeping official plugin adoption separate from third-party marketplace ingestion. Prefer official or clearly maintained plugins, record provenance, review manifests, dependencies, permissions, and MCP surfaces, then compare community marketplace options only after the official path has been evaluated.

- **Official trusted-source lane**
  Use Anthropic's official Claude Code plugins directory as a trusted-source lane for plugin selection: verify provenance against the official repository, inspect install instructions, manifests, dependencies, permissions, and MCP surfaces before installation, check for plugin-skill name, trigger, command, and tool conflicts, and allow official plugins to supersede community marketplace packages when they provide the same capability and pass local review.

- **Anthropic-managed plugin provenance**
  Use the official Anthropic-managed Claude Code plugins directory as the first review source for matching plugin needs before third-party marketplaces. Record the official repository provenance, pin the selected version or commit when supported, inspect dependency and permission changes, check for namespace, command, MCP, and trigger conflicts, and fall back to community marketplaces only when the official directory lacks the required capability or fails the local review criteria.

5. **Agent package-manager evaluation**
   Evaluate agent skill managers and package managers such as AKM by reviewing the trust model, dependency pinning support, manifest schema, cross-agent compatibility, uninstall hygiene, namespace conflict handling, and reproducible install behavior before adopting them for shared skill, command, tool, or knowledge distribution.

- **AKM-style cross-agent kit management**
  Apply package-manager workflows for cross-agent bundles of skills, commands, tools, and knowledge: evaluate manifests and provenance, resolve dependencies, generate and verify lockfiles, pin versions, detect namespace and dependency conflicts, perform integrity checks, preserve reproducible installation, define rollback procedures, and validate portability across supported agent environments.

6. **Trust and provenance review**
   Record source URL, maintainer identity, repository age, release cadence, stars, issues, license, commit activity, and whether installation instructions are traceable to the upstream project.

7. **Dependency and tool-risk inspection**
   Review package files, scripts, MCP declarations, shell commands, network calls, and requested tool permissions before installation. Flag broad filesystem access, credential handling, destructive commands, hidden downloads, or unpinned dependencies.

8. **Namespace and trigger conflict analysis**
   Compare candidate names, descriptions, commands, MCP server IDs, and skill trigger wording against the local skill/plugin registry to prevent ambiguous activation or accidental override of trusted internal workflows.

9. **Installation hygiene**
   Prefer isolated evaluation first. Install into a scratch or staging environment when possible, avoid bulk importing entire catalogs, pin versions or commit SHAs when supported, and keep a rollback path for every adopted package.

10. **Adoption decision record**
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
- Source finding: https://github.com/itlackey/akm
