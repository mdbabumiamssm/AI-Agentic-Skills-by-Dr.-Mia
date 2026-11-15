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
name: 'ai-rulez-agent-configuration-governance'
description: 'Govern AI agent rules, skills, plugins and MCP configuration with ai-rulez through canonical authoring, reproducible generation, artifact validation, dependency locking and provenance checks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# AI Rulez Agent Configuration Governance

## Overview

Use ai-rulez to maintain a canonical source for agent instructions and configuration, generate assistant-specific artifacts, and check whether changes are reproducible and traceable. This skill separates configuration integrity and provenance from behavioral evaluation: a valid or signed configuration does not establish that an agent follows it or performs tasks correctly.

## When to Use This Skill

- Consolidate duplicated assistant rules, skills or MCP settings.
- Detect generated-file drift during review or continuous integration.
- Review changed remote includes, installed skills or plugin dependencies.
- Establish lock and signature checks for agent configuration releases.
- Prepare plugin artifacts for an explicitly requested distribution workflow.

## Core Capabilities

1. **Identify the applicable configuration contract.** Read the existing configuration, selected assistants, repository instructions and deployment scope. Record the installed CLI version and consult `ai-rulez --help` plus command-specific help before using flags. Pin the same release in development and CI; do not assume current website examples match an older binary.

2. **Maintain canonical source content.** For v5 projects, use `.ai-rulez/config.toml` and the associated source directories. Edit canonical inputs rather than patching generated guidance. For a new project, `ai-rulez init acme` initializes a configuration; choose a real project name and only the required presets. Existing projects should retain their chosen profiles and source organization. The documented v5 configuration uses this shape:

   ```toml
   version = "5.0"
   name = "acme"
   presets = ["claude", "cursor", "codex"]
   ```

   Configuration details and initialization are documented in the [project overview](https://goldziher.github.io/ai-rulez/).

3. **Verify target formats and coverage.** Check the [supported harness matrix](https://goldziher.github.io/ai-rulez/harnesses/) for each selected preset. Documented examples include `AGENTS.md`, `CLAUDE.md`, Cursor `.mdc` rules, skill directories, JSON MCP settings and Codex TOML settings. Support differs by feature: hooks, permissions and path-scoped rules do not necessarily retain identical semantics across assistants. Review omitted features, fallback behavior and merged settings explicitly; the advertised assistant count is not a compatibility guarantee.

4. **Generate and validate artifacts.** Run `ai-rulez generate`, inspect the resulting diff, and run `ai-rulez generate --check` to detect drift. Run `ai-rulez validate --strict` and investigate findings using the reported rule codes; `ai-rulez validate --explain AR001` illustrates the explanation interface. Also parse generated JSON or TOML with suitable parsers and inspect Markdown/frontmatter against the destination format. Verify MCP executable names, arguments, environment references and endpoints without starting servers merely to validate configuration. Preserve user-managed keys in merged documents. See [strict validation](https://goldziher.github.io/ai-rulez/strict-validation/).

5. **Lock dependencies and review changes.** Use `ai-rulez lock --diff` for review, then `ai-rulez lock` to accept intended pins in `.ai-rulez/ai-rulez.lock`. Use `ai-rulez lock --check` as the gate; `--diff` exits successfully even when differences exist. Review source revisions and content changes before accepting a refresh. With reviewed pins, use `ai-rulez generate --locked`; `--frozen` additionally avoids network access and requires available dependencies. Keep the release, profile and inputs fixed when checking reproducibility. See the [lock-file guide](https://goldziher.github.io/ai-rulez/lockfile/).

6. **Check provenance separately.** When signing is in scope, configure the documented `[signing]` policy for the intended signer identity and issuer. After final lock and approval changes, use `ai-rulez sign --lock`; verify with `ai-rulez verify --attestation` and retain `ai-rulez lock --check`. A signature over the lock does not independently establish that working sources match it. Record verification results and the exact subject covered; plugin bundles require separate treatment. See [signing the lock](https://goldziher.github.io/ai-rulez/lockfile/#signing-the-lock).

7. **Separate optional packaging and evaluation.** The repository documents `generate --plugin`, `verify --plugin`, `publish`, `review` and `eval`. Resolve their required arguments from the installed version before use. Package and verify locally when requested; publishing requires an authorized destination. For behavioral evaluation, define representative tasks, expected outcomes, model/harness versions and retained run evidence. Report evaluation results separately from lint, drift, lock and signature checks. See the [repository command overview](https://github.com/Goldziher/ai-rulez).

### Execution and completion criteria

Begin with a bounded inventory, then author, generate, inspect, reconcile intended pins and run the final checks. If an enforced lock blocks intentional edits, follow the installed version's documented update procedure; do not disable enforcement or rewrite pins merely to make a gate pass.

Complete the configuration workflow when validation, generation drift checks and the applicable lock/signature gates pass, with no unexplained artifact changes. Record command exit statuses and unsupported target features. On missing dependencies, incompatible versions or failed provenance checks, report the specific blocker and remaining work instead of claiming completion. The 15-minute outcome is a workflow target, not a measured performance guarantee.

## Inputs / Outputs

**Inputs:** Repository root; canonical rules, skills and configuration; selected assistant presets and profiles; pinned CLI release; dependency references and existing lock; applicable organization policy; and signer identity or evaluation requirements when requested.

**Outputs:** Reviewed canonical changes and generated artifacts; a target-format coverage summary; validation and drift results; reviewed dependency pins; and a concise pass/fail report with unresolved findings. Include signature evidence, verified plugin bundles or behavioral evaluation records only when those operations are in scope. Keep credentials out of source files and reports, using the destination's supported secret references.

## References

- [ai-rulez source repository and command overview](https://github.com/Goldziher/ai-rulez)
- [Canonical configuration and v5 overview](https://goldziher.github.io/ai-rulez/)
- [Supported harnesses and generated formats](https://goldziher.github.io/ai-rulez/harnesses/)
- [Strict validation](https://goldziher.github.io/ai-rulez/strict-validation/)
- [Lock files, reproducibility and signing the lock](https://goldziher.github.io/ai-rulez/lockfile/)
