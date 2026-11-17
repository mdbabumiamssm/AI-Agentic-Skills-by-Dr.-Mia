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
name: 'atlas-coding-agent-change-tracking'
description: 'Connect supported coding agents to Atlas, trace sessions to repository changes, query their context, and prepare attributed changes for Git review.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Atlas Coding Agent Change Tracking

## Overview

Use Atlas to connect coding-agent activity with repository changes and recover the session context behind a commit. This skill produces a reviewable attribution record and Git handoff while distinguishing captured evidence from inferred authorship. Atlas offers desktop sessions and a separate terminal integration; select the route appropriate to the user's environment. [Atlas repository](https://github.com/tryatlas/atlas)

## When to Use This Skill

- Connect existing coding agents to a repository's Atlas capture workflow.
- Determine which sessions contributed to a commit and inspect their instructions or tool calls.
- Review changes from several agents before staging or opening a pull request.
- Investigate missing attribution after imports, rebases, amends, or squash merges.
- Prepare a handoff that preserves session context alongside the team's normal Git review.

## Core Capabilities

1. **Select a supported integration.** Desktop Atlas prioritizes Claude Code and Codex through ACP and offers additional registry agents, including OpenCode and Kilo Code; upstream notes that testing of other registry agents is ongoing. The CLI documents Claude Code, Codex, Cursor, Copilot, OpenCode, Droid, Antigravity, and Pi. Do not infer CLI support from desktop registry availability. [Desktop integrations](https://github.com/tryatlas/atlas), [CLI integrations](https://docs.tryatlas.cc/docs/cli/agents)

2. **Establish repository scope and capture.** Record the checkout, branch, starting commit, and pre-existing modifications before enabling capture. Keep unrelated user changes distinguishable from the task. Choose local storage unless sharing is requested; cloud setup must use the intended organization and project.

   Inspect the baseline from the intended repository:

   ```bash
   git rev-parse --show-toplevel
   git branch --show-current
   git rev-parse HEAD
   git status --short
   ```

   For an installed CLI, inspect its version and current binding first:

   ```bash
   atlas-cli --version
   atlas-cli status
   atlas-cli init --help
   ```

   For a new local setup, use `atlas-cli init --local`. For requested cloud setup, use `atlas-cli login` followed by `atlas-cli init` and select the correct organization, project, agents, and sync setting. Do not reinitialize a working binding merely to inspect it. Verify version-specific flags with command help. [Command reference](https://docs.tryatlas.cc/docs/cli/commands)

   Desktop alternative: open the repository, open Session capture from the project breadcrumb, select Local, and enable it. Cloud capture adds sharing to the local record; inspect the history-sharing preview before publishing. Capture is initially disabled per project. [Session capture](https://docs.tryatlas.cc/docs/source-control/capture)

3. **Connect agents and verify a captured task.** On the CLI route, add only the requested integrations, for example:

   ```bash
   atlas-cli agents add claude-code
   atlas-cli agents add codex
   atlas-cli agents list
   atlas-cli capture enable
   atlas-cli agents sessions
   atlas-cli status
   ```

   These integrations modify project agent configuration to register hooks and shared memory. Inspect those changes and complete the agent's own trust prompts as appropriate. Launch the selected agent normally, perform one small authorized task, then repeat the session and status checks. If installation is missing, use the upstream CLI installation instructions; do not substitute an unrelated package named Atlas. [CLI setup](https://github.com/tryatlas/cli)

   Desktop alternative: create a session, select its agent in the composer, and give it the scoped task. Record the session ID, project, and branch; keep separate tasks identifiable. Additional agents are optional, not a requirement to launch parallel work. [Sessions](https://docs.tryatlas.cc/docs/product/chat)

4. **Review and query attribution.** Locate the task's session, inspect its prompts and tool calls, and compare the affected paths against the actual Git diff. After a commit, open Source Control history and follow its session links. Desktop checkpoint chat can answer questions about a selected recorded session; treat generated explanations as leads to verify against the transcript and patch. [Git review](https://docs.tryatlas.cc/docs/source-control/git), [Checkpoint querying](https://github.com/tryatlas/atlas)

   The CLI also exposes shared-memory browsing with `atlas-cli memory browse`; connected agents can use the registered memory tools, including `memory_search` and `memory_why`, through their discovered tool schemas. These are MCP tools, not invented shell subcommands. Record the evidence supporting each answer and mark unresolved attribution explicitly. [CLI memory](https://github.com/tryatlas/cli)

   Attribution has limits: imported sessions lack checkpoints; edits to existing files may remain linked after human rewrites. Desktop capture observes Git history, whereas CLI capture uses a post-commit hook. History rewriting may orphan links when no unambiguous matching patch remains. A checkpoint is evidence of participation, not proof of exclusive line authorship. [Capture semantics](https://docs.tryatlas.cc/docs/source-control/capture)

5. **Integrate a reviewed change with Git.** Inspect both working-tree and staged changes, including untracked files shown by status:

   ```bash
   git diff
   git diff --cached
   git diff --check
   ```

   Stage only the intended files or hunks, run the repository's relevant checks, and review the staged patch again. Follow the team's existing branch and commit conventions; commit, push, or publish a pull request only within the requested scope. Include session references and validation results in the review handoff. Atlas's turn diff excludes edits made through shell commands, so use Git's complete patch for final review. [Diff scopes](https://docs.tryatlas.cc/docs/source-control/git)

   Once a commit is authorized and created, record its SHA and confirm the session link in Atlas. Recheck attribution after any history rewrite; preserve missing or orphaned status rather than manufacturing a replacement association.

6. **Diagnose incomplete capture.** Inspect `atlas-cli status` and the relevant command's `--help`, or the desktop capture indicator's reason. Check project binding, agent hooks, and storage permissions before repeating the small verification task. Desktop capture can be degraded by another recording window or an unwritable store. Stop repeated retries when the same cause persists and report the evidence needed to resolve it. [CLI diagnostics](https://docs.tryatlas.cc/docs/cli/commands), [Capture troubleshooting](https://docs.tryatlas.cc/docs/source-control/capture)

## Inputs / Outputs

**Inputs:** Repository path; requested task or existing commit to investigate; selected agents; installed Atlas version and route; branch/base commit; local or cloud preference; existing project binding; and the team's review and validation commands.

**Outputs:** A concise record containing:

- Repository, branch, baseline SHA, Atlas version, storage mode, and connected agents.
- Session identifiers, affected paths, checkpoint/commit associations, and unresolved attribution.
- Reviewed patch scope, check results, and outstanding conflicts or capture failures.
- Commit SHA or pull-request reference when created within scope; otherwise a ready-to-review handoff.

**Completion check:** Within the 15-minute workflow target, verify at least one relevant recorded session, compare its evidence with the repository patch, and produce the handoff. A commit-linking task additionally requires an observed checkpoint association. If installation, authentication, capture, or attribution blocks verification, report partial completion and the precise blocker instead of claiming success; the time target is not a performance benchmark.

## References

- [Atlas repository and desktop overview](https://github.com/tryatlas/atlas)
- [Atlas CLI repository, installation, and memory tools](https://github.com/tryatlas/cli)
- [CLI supported agents](https://docs.tryatlas.cc/docs/cli/agents)
- [CLI command reference](https://docs.tryatlas.cc/docs/cli/commands)
- [Sessions](https://docs.tryatlas.cc/docs/product/chat)
- [Session capture and checkpoint limitations](https://docs.tryatlas.cc/docs/source-control/capture)
- [Git and diff review](https://docs.tryatlas.cc/docs/source-control/git)
