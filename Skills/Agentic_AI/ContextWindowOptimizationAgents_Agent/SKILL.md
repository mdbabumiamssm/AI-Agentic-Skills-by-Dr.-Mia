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
name: 'context-window-optimization-agents'
description: 'Optimize AI coding-agent context by sandboxing verbose tool output, preserving useful summaries, and managing transcript hygiene across agent runtimes.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Context Window Optimization for Agents

**Overview**

This skill helps AI coding and research agents reduce context-window waste from verbose tool output while preserving the evidence needed to make correct decisions. It applies context-mode-style sandboxing, summary extraction, transcript hygiene, and cross-runtime context budgeting for agent workflows that involve logs, searches, builds, tests, and repository exploration.

Use it to keep the active context focused on task state, decisions, errors, changed files, and next actions instead of raw command output.

**When to Use This Skill**

- Tool outputs are long enough to crowd out the user request, plan, code context, or current findings.
- Build, test, lint, search, or agent-runtime logs need to be summarized without losing actionable failures.
- A coding agent must preserve useful context across multiple commands, files, or platforms.
- The workflow involves context-mode, Codex, Claude Code, Cursor, Copilot, OpenCode, Zed, MCP tools, or similar agent runtimes.
- A transcript needs cleanup before handoff, compaction, replay, or escalation to another agent.
- The user asks for context budgeting, token reduction, output sandboxing, agent memory hygiene, or tool-output compression.

**Core Capabilities**

1. **Context budget assessment**: Estimate which parts of the transcript and tool output are still decision-relevant, stale, duplicated, or safe to summarize.
2. **Tool-output sandboxing**: Keep raw verbose output outside the active reasoning path when possible, and surface only the command, status, key lines, errors, file paths, and implications.
3. **Actionable summarization**: Convert large logs into concise findings that retain failure modes, stack traces, test names, changed files, and reproduction commands.
4. **Transcript hygiene**: Remove repeated discoveries, outdated plans, superseded assumptions, and low-value narration from handoff summaries.
5. **Cross-agent handoff packaging**: Produce compact state bundles that include objective, constraints, current files, decisions made, open risks, and next steps.
6. **Platform-aware adaptation**: Adjust the workflow for agent runtimes that support hooks, plugins, MCP servers, CLI wrappers, or local transcript storage.
7. **Verification preservation**: Keep enough raw anchors, command names, exit states, and file references for another agent or human to audit the summary.
8. **Context-mode pattern checks**: Sandbox verbose tool output, preserve compressed summaries with decision-critical evidence, verify coverage across supported agent platforms, and measure context-token reduction without dropping actionable signals.
9. **Concrete implementation reference**: Use `mksglu/context-mode` as a TypeScript reference for coding-agent transcript hygiene: sandbox tool output, keep compact summaries active, adapt hooks/plugins/MCP integrations across supported platforms, document failure modes such as omitted error evidence or stale summaries, and preserve raw logs when stack traces, failing commands, reproduction details, or audit requirements matter more than compression.
10. **Raw-log promotion rules**: Measure token reduction after sandboxing verbose output, keep compact summaries in active context, and promote raw logs back into working context when summaries omit error evidence, stack traces, failing commands, reproduction details, or other decision-critical signals.
11. **Sandbox replay pointers**: Isolate verbose logs behind replay pointers while preserving compact summaries, compare observed token savings with the source-reported 98% reduction claim, and recover raw output when debugging requires full command output, stack traces, failing test logs, or audit evidence.
12. **Cross-runtime context-mode integration**: Apply context-mode patterns for Codex, Claude Code, Cursor, Copilot, and MCP servers by sandboxing verbose tool output, keeping concise summaries and artifacts in active context, measuring token reduction, and defining hidden-output failure modes such as missing stack traces, failed command context, stale summaries, or unavailable raw-log replay.
13. **High-signal context-mode hygiene pattern**: Treat `mksglu/context-mode` as a coding-agent pattern for sandboxing tool output, compacting transcripts, supporting many agent platforms, tracking the source-reported 98% token reduction, and explicitly guarding against hidden diagnostically important output.
14. **Context-mode operational controls**: For coding-agent runs, sandbox verbose tool output, keep compact summaries active, measure token reduction before and after filtering, define raw-log recovery paths for omitted diagnostics, and document cross-runtime compatibility checks for the agent platforms in use.
15. **Context-output sandbox evaluation**: Compress tool output into durable summaries for cross-agent reuse, then compare token reduction against lost debugging signal by checking whether errors, stack traces, failing commands, reproduction details, stale summaries, and raw-output recovery paths remain available.
16. **High-adoption context-mode example**: Use `mksglu/context-mode`—reported in the source finding with 12,504 GitHub stars and support for 14 platforms—as an example of sandboxing verbose coding-agent tool output; verify summary fidelity and transcript hygiene, test compatibility with the target platform, document filtering and replay failure modes, and disable suppression when it could hide diagnostic evidence needed for debugging or audit.
17. **Installation-to-evaluation operating procedure**: For each supported coding agent, use its documented native plugin path when available or register Context Mode as an MCP server with the platform's hooks and routing instructions; restart, verify connectivity and hook activation, and avoid duplicate plugin, MCP, or hook registrations. Route verbose shell, file, search, web, test, and build output through sandbox tools while retaining compact summaries and raw-output recovery pointers. Evaluate the same representative task before and after by recording raw versus retained token or byte counts, checking task correctness and summary fidelity, and confirming that errors, stack traces, failing commands, user decisions, and compaction recovery remain available; treat the source-reported 98% reduction as a claim to reproduce, not a universal result. Fall back to raw output when hooks fail, routing is unavailable, summaries become stale, or diagnostics disappear. Confirm local storage, retention and purge behavior, secret handling, inherited permission rules, and network/data-access boundaries before processing private material.
18. **Coding-agent context hygiene adoption example**: Treat `mksglu/context-mode` as a high-adoption TypeScript example for sandboxing verbose tool output, retaining structured summaries, measuring token reduction, and supporting multiple agent runtimes; define recoverability requirements for suppressed raw output when summaries omit errors, stack traces, failed commands, reproduction details, stale-summary conflicts, or audit evidence.
19. **Auditable context-mode compaction pattern**: Apply the `mksglu/context-mode` finding as a coding-agent transcript hygiene pattern: sandbox verbose tool output, keep compact summaries plus command, status, error, and file anchors in active context, integrate through platform-native plugins, hooks, or MCP where available, measure retained versus raw token or byte counts per run, and fail open to raw output when compaction hides stack traces, failed-command context, reproduction details, stale-summary conflicts, or audit evidence.
20. **Context-mode implementation adoption check**: Use `mksglu/context-mode` as a concrete implementation pattern for sandboxing verbose tool output across coding-agent runtimes; preserve transcript hygiene with compact summaries and raw-output recovery anchors, verify the source-reported 14-platform coverage against the target runtime, and complete integration-risk checks for hooks/plugins/MCP registration, storage and retention, secret handling, permissions, and failure-open behavior before adoption.
21. **Audit-ready transcript hygiene budgets**: For coding-agent context-mode workflows, sandbox verbose tool output, summarize retained evidence, preserve command provenance, set transcript hygiene budgets for what stays active versus recoverable, and validate that reduced context still supports auditability and recovery.
22. **Context-mode coding-agent hygiene pattern**: Use `mksglu/context-mode` as a concrete TypeScript implementation pattern for sandboxing verbose tool output, preserving concise summaries, defining hooks/plugins/MCP or platform-native integration points across agent runtimes, measuring observed token reduction against the source-reported 98% reduction claim, and documenting failure modes where hidden output removes needed evidence such as errors, stack traces, failing commands, reproduction details, stale-summary conflicts, or audit anchors.
23. **Cross-agent context-mode adoption pattern**: Treat `mksglu/context-mode` as a concrete cross-agent optimization pattern: route verbose tool output into sandboxed, recoverable storage; keep only concise summaries, command status, file paths, and decision-critical evidence in the active transcript; verify the target runtime against the source-reported 14-platform support before adoption; regression-test hidden-output behavior for missing errors, stale summaries, or lost reproduction context; and preserve raw logs whenever stack traces, failing commands, full test/build output, audit evidence, or summary-fidelity checks are required.
24. **Concrete context-mode transcript-hygiene implementation**: Treat `mksglu/context-mode` as a TypeScript implementation pattern for coding-agent transcript hygiene: sandbox verbose tool output, retain compact summaries plus command, status, and file anchors, verify the source-reported 14-platform coverage for the target runtime, preserve raw logs when stack traces, failing commands, reproduction details, full test or build output, audit evidence, or summary-fidelity checks are needed, and evaluate the source-reported 98% context reduction against any added debugging cost.
25. **Context Mode token-reduction validation**: Use `mksglu/context-mode` as a concrete coding-agent pattern for sandboxing verbose tool output while retaining summaries, command status, file anchors, and recoverable raw logs; verify cross-runtime support against the source-reported 14 platforms, measure observed token reduction instead of assuming the reported 98% result, and fail open when summaries hide errors, stale state, stack traces, failing commands, reproduction details, or audit evidence.
26. **High-signal sandbox implementation reference**: Use `mksglu/context-mode` as an implementation reference for sandboxing verbose coding-agent tool output, preserving compact summaries, reducing context usage across coding-agent platforms, and keeping transcripts hygienic while retaining recoverable raw evidence for debugging.
27. **Concrete cross-agent context optimization pattern**: Use `mksglu/context-mode` across Codex, Claude Code, Cursor, Copilot, and MCP tool workflows by sandboxing verbose tool output, retaining compact summaries with command status and file anchors, preserving auditability through recoverable raw-output pointers, measuring observed token reduction against the source-reported 98% claim, and documenting per-runtime integration considerations for hooks, plugins, MCP registration, permissions, and failure-open behavior.
28. **Context-mode transcript hygiene guardrails**: Sandbox verbose tool output outside the active transcript, preserve compact summaries with command status and evidence anchors, quantify raw versus retained token or byte counts against the source-reported 98% reduction claim, confirm compatibility constraints for each target agent runtime in the source-reported 14-platform scope, and keep recoverable raw evidence when auditability, stack traces, failing commands, or reproduction details could be lost.
29. **Cross-platform verbose-output sandbox pattern**: Apply `mksglu/context-mode` as a concrete pattern for Codex, Claude Code, Cursor, Copilot, and MCP workflows: sandbox verbose tool output, preserve compact active summaries, measure observed token reduction against raw output, and maintain transcript hygiene without assuming the source-reported 98% reduction will reproduce in every run.
30. **High-signal implementation reference**: Treat `mksglu/context-mode` as a TypeScript reference for coding-agent transcript hygiene, sandboxed tool output, aggressive summarization, source-reported 14-platform support, and measurable context reduction while retaining actionable debugging evidence.
31. **High-adoption context hygiene reference**: Evaluate `mksglu/context-mode` as a high-adoption implementation reference for coding-agent context hygiene: sandbox verbose tool output, retain concise summaries in active context, document which source-reported supported platforms apply to the target workflow, and validate measurable token reduction before adoption instead of assuming the source-reported 98% result.
32. **Context-output sandboxing contract**: Use `mksglu/context-mode` as the concrete context-output sandboxing pattern for the source-reported 14-platform coverage; validate the reported 98% reduction on representative local tasks, require summaries to retain command, exit status, file paths, failing test names, error messages, stack-trace anchors, reproduction commands, and raw-output recovery pointers, preserve raw logs for full stack traces, flaky failures, audit evidence, privacy review, and summary-fidelity checks, and add regression tests that prove sandboxing does not lose diagnostic detail needed to reproduce or fix failures.

**Inputs / Outputs**

Inputs:

- User objective, current task constraints, and any platform/runtime requirements.
- Raw or summarized tool output from shell commands, code search, tests, linters, builds, web fetches, MCP tools, or agent hooks.
- Repository state such as relevant files, changed paths, test commands, errors, and prior decisions.
- Optional context budget target, such as a maximum token budget, summary length, or handoff format.

Outputs:

- A compact context summary with current objective, important facts, command results, failures, changed files, and next actions.
- A classification of raw output into keep, summarize, archive, or discard categories.
- A context-safe handoff note for another agent, future session, or transcript compaction.
- Optional runtime guidance for configuring tool-output filtering, command wrappers, hooks, or MCP-based context management.

**References**

- Source finding: https://github.com/mksglu/context-mode
- Model Context Protocol specification: https://github.com/modelcontextprotocol/specification
