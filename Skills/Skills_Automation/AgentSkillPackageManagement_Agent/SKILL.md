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
name: 'agent-skill-package-management'
description: 'Install and manage versioned AI agent skills, commands, tools, and knowledge bundles reproducibly with manifests, trust checks, lockfiles, updates, rollback, and safe removal.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Agent Skill Package Management

## Overview

Use this skill to install and manage agent capability packages across supported AI coding agents and runtimes. Apply reproducible manifests, pinned versions, provenance checks, dependency controls, lockfiles, and reversible lifecycle operations so package state remains auditable and consistent.

## When to Use This Skill

- Install a skill, command, tool, MCP integration, or knowledge bundle from a registry or repository.
- Create or validate a package manifest before publishing or deployment.
- Pin package and dependency versions for reproducible environments.
- Verify package origin, integrity, licenses, signatures, checksums, or repository trust signals.
- Resolve dependency constraints and detect incompatible or duplicated packages.
- Deploy the same package set across multiple agent platforms or workspaces.
- Generate, refresh, compare, or enforce a lockfile.
- Preview and apply package updates without silently changing behavior.
- Roll back a failed installation or update to a known-good state.
- Remove a package without leaving orphaned dependencies, files, or configuration.

## Core Capabilities

1. **Inventory the target environment.** Identify the package manager, agent runtimes, configuration roots, installed packages, active versions, local modifications, and platform-specific deployment paths before making changes.

2. **Validate manifests.** Require a machine-readable manifest containing package identity, version, package type, source, supported agents, entry points, dependencies, permissions, files, and integrity metadata. Reject ambiguous names, malformed versions, unsupported targets, and undeclared executable behavior.

3. **Pin versions and sources.** Prefer immutable releases, commit hashes, content digests, or exact semantic versions. Record the canonical source URL and resolved revision; avoid floating branches or unbounded version ranges in reproducible deployments.

4. **Check provenance and trust.** Confirm that the package resolves to the expected repository or registry, inspect ownership and release history, verify available signatures or checksums, review requested permissions, and flag install scripts or binaries that require additional scrutiny. Do not treat popularity as proof of safety.

5. **Resolve dependencies.** Build the complete dependency graph, apply version constraints, detect cycles, separate required and optional dependencies, and fail with a clear explanation when no compatible solution exists. Do not silently downgrade or substitute packages.

6. **Apply namespacing.** Use stable package identifiers and target-aware namespaces for commands, tools, skills, and knowledge resources. Preserve aliases only when they are explicit and do not shadow existing capabilities.

7. **Detect conflicts.** Check command names, tool identifiers, file destinations, configuration keys, environment variables, ports, dependency versions, and agent-specific registration entries. Present conflicts and the proposed resolution before changing state.

8. **Plan cross-agent deployment.** Translate the package into each target agent's supported directory and configuration format while preserving a shared package identity and pinned source. Report unsupported features instead of pretending deployments are equivalent.

9. **Create and enforce lockfiles.** Record direct and transitive package versions, source revisions, integrity digests, target platforms, and relevant resolver metadata. Treat the lockfile as the authoritative resolved state and update it atomically with package changes.

10. **Install transactionally.** Produce a dry-run plan, back up affected state, stage files in a temporary location, validate them, and then commit the installation. On failure, restore the prior state and retain actionable diagnostics.

11. **Update deliberately.** Compare installed, locked, and available versions; review manifest, dependency, permission, and file changes; then update selected packages under explicit constraints. Regenerate the lockfile and verify all target agents after the update.

12. **Roll back safely.** Restore the previous package files, manifest resolution, lockfile, configuration, and registrations from a recorded transaction or known-good revision. Verify that rollback removes partial state introduced by the failed change.

13. **Remove with dependency awareness.** Identify reverse dependencies, block unsafe removal, unregister agent integrations, delete only package-owned files, and prune dependencies only when no remaining package requires them. Preserve user-authored or locally modified content unless removal is explicitly approved.

14. **Verify and report.** Confirm package discovery, entry-point loading, dependency integrity, target-agent registration, and lockfile consistency. Return a concise report of resolved versions, provenance, changes, warnings, rollback data, and verification results.

15. **Compare AKM package-manager patterns.** Use AKM (Agent Kit Manager) as another package-manager example for AI agent skills, commands, tools, and knowledge bundles. Compare AKM with ccpi/craftdesk-style managers by checking manifest structure, reproducible installs, pinned source and version handling, source trust and provenance review, lock/update behavior, and cross-agent portability.

16. **Compare Craftdesk package-manager patterns.** Include Craftdesk as a TypeScript package-manager example for Claude Code skills, agents, and related AI resources. In the comparison matrix, review dependency handling and metadata, install/update/remove hygiene, version tracking, lockfile expectations, rollback expectations, trust review for each third-party bundle or source before adoption, and conflict handling when multiple skill sources target the same skill, agent, command, file, or registration.

## Inputs / Outputs

### Inputs

- Requested operation: inspect, install, deploy, update, roll back, or remove.
- Package identifier, manifest, registry coordinate, repository URL, release tag, or commit.
- Target agent platforms, workspaces, profiles, and destination paths.
- Version constraints, trust policy, permission policy, and allowed package sources.
- Existing package inventory, configuration, and lockfile when present.
- Conflict policy, such as fail, namespace, replace, or request explicit approval.

### Outputs

- Validated package manifest and normalized package identity.
- Resolved dependency graph with exact versions and immutable source references.
- Provenance, integrity, permission, and conflict-check results.
- Dry-run change plan describing files, configuration, registrations, and dependencies.
- Updated package state and lockfile, or an unchanged state when validation fails.
- Cross-agent deployment matrix showing success, unsupported features, and deviations.
- Rollback record or removal report with retained files and dependency decisions.
- Final verification summary with warnings and remediation steps.

## References

- [akm - Agent Kit Manager](https://github.com/itlackey/akm)
- [Craftdesk](https://github.com/mensfeld/craftdesk)
