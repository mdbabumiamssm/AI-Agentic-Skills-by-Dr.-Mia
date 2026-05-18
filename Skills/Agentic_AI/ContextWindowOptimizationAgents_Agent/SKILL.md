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
