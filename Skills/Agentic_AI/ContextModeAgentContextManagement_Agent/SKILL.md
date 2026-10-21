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
name: 'context-mode-agent-context-management'
description: 'Configure and validate Context Mode MCP tools and client hooks to contain tool output, retrieve indexed evidence, and preserve coding-agent session state while measuring context savings and task fidelity.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Context Mode Agent Context Management

## Overview

Use Context Mode to process bulky tool results outside the conversation, retrieve relevant evidence from an index, and recover working state through supported client hooks. This skill guides client-specific setup and a small validation exercise that measures output reduction without losing facts needed for the task. The repository advertises 98% reduction and support across multiple platforms; treat savings and feature coverage as claims to verify against the installed version and workload.

## When to Use This Skill

- Logs, test output, repository searches, web pages, or API responses consume excessive conversation context.
- A coding agent loses decisions, unresolved errors, or task state after compaction or resume.
- Context Mode MCP tools connect, but large outputs still bypass routing.
- A team needs to install, troubleshoot, or compare Context Mode behavior across coding clients.
- Context savings must be demonstrated alongside correct answers and recoverable evidence.

## Core Capabilities

1. **Select a supported integration.** Identify the client, client version, Context Mode version, runtime, configuration scope, and available tools. Read the current repository installation section and platform support matrix; MCP connectivity, interception, event capture, and compaction recovery are distinct capabilities. Do not infer complete support from the repository's platform count.

2. **Install and verify the requested client.** Apply only the integration the user requested, merging configuration with existing entries and preserving unrelated hooks and instructions. Common installation paths from the upstream README are:

   - **Claude Code:** Run `/plugin marketplace add mksglu/context-mode` and `/plugin install context-mode@context-mode` in the client, reload plugins or restart, then run `/context-mode:ctx-doctor`. These are client commands, not shell commands. The MCP-only alternative lacks plugin hook enforcement.
   - **Codex CLI:** Follow the current Codex subsection for marketplace or manual setup. The documented plugin path starts with `codex plugin marketplace add mksglu/context-mode`, then installation through the plugin UI. Where required by the installed client, enable `hooks` and `plugin_hooks` under the existing `[features]` table and complete the client's hook trust flow. Avoid duplicating manual hooks when plugin hooks already provide them.
   - **Other clients:** Use their matching configuration under upstream `configs/` and their current installation instructions. Generic MCP registration does not install lifecycle hooks. For instruction-only clients, merge the appropriate routing guidance into the existing instruction file; do not overwrite project rules.

   Confirm runtime prerequisites from the chosen installation path; several current paths require Node.js 22.5 or newer, or Bun. Restart the client and inspect tool discovery. Run `ctx_doctor` and `ctx_stats` through the exposed MCP tools; the CLI diagnostic equivalent is `context-mode doctor`. Tool schemas and client prefixes govern invocation, so inspect them instead of guessing parameters. The frontmatter tools support setup and inspection; actual MCP validation also requires the host to expose Context Mode tools. If they are absent, report that dependency rather than claiming a completed tool test.

3. **Route bulky work before it enters context.** Use `ctx_execute` for code that computes a compact answer, `ctx_execute_file` for analysis of a local file, and `ctx_batch_execute` for related commands and queries. Print the requested findings, error details, counts, and source locations rather than entire inputs. Keep ordinary edits in the host's appropriate editing tools. Each execution has separate process state; do not depend on variables surviving another call. Output containment is not a full operating-system security boundary, and execution must respect existing filesystem and command permissions.

4. **Index and retrieve evidence progressively.** Use `ctx_index` for content to revisit and `ctx_fetch_and_index` for URL content, then `ctx_search` for focused queries. Preserve returned source identifiers and search terms so findings can be traced to indexed content. Ask for missing details using targeted follow-up queries; an empty or narrow result does not prove the source lacks the information. If freshness matters, use the installed tool's documented cache bypass, such as `force: true` or `ttl: 0`, and record the refresh. Do not paste a full document into the conversation merely to index it afterward.

5. **Validate routing independently of connectivity.** Confirm that the client loads the configured hooks and use a harmless, output-producing read or command covered by the routing rules. Observe whether the call is blocked, redirected, or simply guided by instructions. Record that observed behavior separately from MCP health. In the upstream Codex integration, blocking support does not imply tool-input rewriting; verify the actual client capability before expecting automatic rerouting. If routing fails, inspect hook registration, matchers, executable paths, and trust state before changing configuration.

6. **Manage session memory deliberately.** Separate the searchable content index from session events and recovery snapshots. On a client with the necessary capture and restore hooks, record a harmless task decision and unfinished step, then perform a supported compaction or resume and verify both are recovered correctly. Use the client's documented resume mechanism; do not assume a fresh session retains the prior session's state. Check actual storage and retention behavior for the installed adapter, and reconcile recovered state with current files and the user's latest instructions. Treat `ctx_purge` as permanent deletion of indexed content, not a universal reset of all session memory; use it only when the requested cleanup includes that scope.

7. **Measure savings and fidelity within a bounded trial.** Aim to complete discovery, setup checks, one representative output test, and a recovery check within 15 minutes. Select a harmless existing log or test result with known facts, including an error or uncommon detail. Determine the baseline answer and raw output size through local processing, without injecting the full baseline into the conversation. Run the equivalent Context Mode analysis and retrieve enough evidence to verify the answer.

   Capture before-and-after `ctx_stats` values and compare the same workload. Report raw bytes, total returned bytes including follow-up retrievals, tool calls, elapsed time, and whether every required fact was recovered. Compute output reduction as `100 * (1 - returned_bytes / raw_bytes)` only when `raw_bytes > 0` and both measurements use the same encoding and scope. Label token estimates separately; byte reduction is not total conversation-token reduction. Include routing and hook overhead when measuring total context cost. Do not substitute the repository's advertised percentage for an observation.

   Success requires reachable tools, an evidenced routing result, a correct compact answer with traceable supporting content, and a measured reduction or an honest finding of no reduction. Report recovery as passed, failed, unsupported, or not tested. If setup or recovery cannot finish within the trial, return a partial result with the exact blocker and next action; do not mark the workflow successful or upgrade software merely to satisfy the time target.

## Inputs / Outputs

**Inputs:** Target client and version; requested installation or troubleshooting scope; existing MCP, hook, and routing configuration; accessible runtime; a representative non-sensitive source; expected facts or task acceptance criteria; and any session recovery or retention requirements.

**Outputs:** A verified configuration or actionable configuration change; diagnostic status; observed routing behavior; compact findings linked to source evidence; session recovery status; and a measurement record containing workload, units, raw and returned sizes, reduction calculation, retrieval overhead, elapsed time, and fidelity result. Distinguish observed results from upstream claims and list any untested capabilities. Do not claim that all platforms provide equivalent memory or routing behavior.

## References

- [Context Mode source repository](https://github.com/mksglu/context-mode) — source finding and implementation.
- [Installation, tools, memory lifecycle, and diagnostics](https://github.com/mksglu/context-mode/blob/main/README.md) — current operational instructions; verify against the installed release.
- [Platform support matrix](https://github.com/mksglu/context-mode/blob/main/docs/platform-support.md) — client-specific integration and hook capabilities.
- [Upstream benchmark report](https://github.com/mksglu/context-mode/blob/main/BENCHMARK.md) — repository-reported measurements, not guaranteed outcomes for another workload.
