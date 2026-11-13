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
name: 'grimoire-oci-agent-config-packaging'
description: 'Package and distribute AI-agent configuration with Grimoire OCI artifacts, digest-pinned installs, reviewed updates, rollback, and runtime compatibility checks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Grimoire OCI Agent Configuration Packaging

## Overview

Use Grimoire (`grim`) to package and distribute skills, rules, agents, MCP server configuration, and bundles through OCI registries with digest-pinned lockfiles. This workflow makes configuration delivery reproducible while checking each runtime's supported formats and preserving a recoverable installation. Content integrity and trust in package behavior require separate evidence.

## When to Use This Skill

- Publish reusable agent configuration to an approved OCI registry.
- Reproduce a reviewed configuration on another workstation or CI runner.
- Review upstream changes before updating locked packages.
- Recover a previous configuration after an unsuccessful update.
- Check whether a package works across selected coding-agent runtimes.

## Core Capabilities

1. **Establish scope and compatibility.** Record the installed `grim` version, artifact kinds, selected runtimes and versions, project or global scope, and registry destination. Check `grim --help` and relevant subcommand help against the documentation for that version before using flags. Consult the [compatibility matrix](https://grimoire.rs/clients.html) for each runtime/artifact pair; support for one kind does not imply support for every kind. Record unsupported or skipped entries explicitly, and treat a missing required capability as a failed acceptance check.

2. **Build a reviewable package.** Inspect the source and included supporting files before building. Exclude credentials, local caches, and unrelated material; examine bundled scripts and MCP launch commands as executable dependencies. Use `grim build <source>` to validate and pack locally. Specify `--kind agent` for agent Markdown and `--kind mcp` for MCP TOML, since their extensions overlap with rules and bundles. Preview `grim release <source> <registry-reference> --dry-run`; inspect all tags it would change. Full semantic versions can also move floating tags; use the documented `--no-cascade` option when only the exact tag should move. Publish with `grim release` only within the requested registry scope, then record the actual registry digest. See the [publishing guide](https://grimoire.rs/publishing.html).

3. **Install the reviewed content.** Before changing a working installation, preserve the declaration, lockfile, installed configuration, and relevant install state. Configure the intended clients and scope before adding anything: `grim add <reference>` declares, locks, and installs. Use `grim init` only for a new project, and `grim install` to materialize an existing reviewed lock. Keep `grimoire.toml` and `grimoire.lock` together in version control. A cold cache still requires registry access; a lockfile is not an offline artifact archive. CI should fail on missing or stale locks rather than silently regenerate them. See the [team and CI guide](https://grimoire.rs/guides/team-ci.html).

4. **Verify digests without conflating integrity and trust.** Record the complete locked OCI reference and compare it with the independently approved manifest digest. Use tooling that verifies the manifest bytes and referenced blobs against their OCI descriptors; stop on any mismatch. Do not compare an extracted directory hash directly with an OCI manifest digest: they represent different objects. Separately review provenance, instructions, scripts, destinations, and requested permissions. A matching digest establishes content identity, not publisher authenticity, harmlessness, or fitness for use. An MCP descriptor pin does not automatically pin a downloaded executable, container tag, package-manager dependency, or remote service; record and constrain those dependencies separately.

5. **Review updates before adoption.** Run `grim status --check --format json` to inspect upstream movement. Its successful exit alone is not a health gate. In an isolated copy with matching configuration, run `grim update`, then compare old and new locks, package contents, rendered files, bundle membership, and MCP commands or endpoints. Check the lock even after failure: an update refused because of local edits can already have advanced its pins. Preserve those edits and resolve the intended version before retrying; do not use `--force` as a routine repair. It can authorize overwrites and update deletions. See the [lifecycle guide](https://grimoire.rs/guides/lifecycle.html).

6. **Rehearse rollback and runtime acceptance.** Restore the known-good declaration and lock together in the isolated environment, then install from those pins without re-resolving tags. Confirm the old digests remain retrievable. Reconcile managed files and install state with the saved inventory; installation alone may leave files introduced by the newer release. Preserve unrelated user files and local edits. For each required runtime, verify configuration parsing, discovery of the expected artifact, and a harmless functional check; test MCP connectivity only within the authorized scope. Report rollback as unverified if it was not exercised.

### Completion criteria

Work on the requested operation; publishing is not required for an installation-only task. Within the 15-minute target, produce a package or installation result plus evidence identifying exact digests, changed paths, compatibility outcomes, and rollback readiness. The target is an execution objective, not a measured benchmark or a guarantee.

Mark success only when required artifacts are present in the intended scope, integrity checks pass, and required runtime checks succeed. If credentials, registry availability, unsupported formats, or unresolved local changes prevent completion, report the specific blocker and completed checks; do not label an untested or partial deployment successful.

## Inputs / Outputs

**Inputs**

- Local artifact sources or OCI references, including bundle members where applicable.
- Existing declaration, lockfile, configuration inventory, and known-good recovery point.
- Requested operation, target runtimes and scope, pinned CLI version, and registry access.
- Approved package digests or provenance policy, intended release tags, and acceptance criteria.

**Outputs**

- For packaging: validated artifact, reviewed publication plan, and published reference/digest if publication was requested and completed.
- For installation or update: paired declaration and lockfile plus the intended runtime configuration changes.
- A concise evidence record containing CLI version, old/new digests, content and configuration changes, integrity results, runtime pass/fail/skipped results, and rollback result or limitation.
- Explicit blockers and remaining checks when the requested operation cannot be completed.

## References

- [Grimoire source repository and project overview](https://github.com/grimoire-rs/grimoire).
- [Publishing skills and rules, artifact kinds, and release previews](https://grimoire.rs/publishing.html).
- [Install/update lifecycle and local edit handling](https://grimoire.rs/guides/lifecycle.html).
- [Digest-pinned team installations and CI](https://grimoire.rs/guides/team-ci.html).
- [Client compatibility matrix](https://grimoire.rs/clients.html).
