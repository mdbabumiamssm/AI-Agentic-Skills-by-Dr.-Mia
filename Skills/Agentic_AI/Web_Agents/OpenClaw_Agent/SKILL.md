---
name: openclaw-agent-operations
description: Operate OpenClaw as a local AI assistant and agent OS with explicit security boundaries. Use when local system automation, browser work, and provider-pluggable execution need to stay on the operator's machine rather than moving to a hosted browser agent surface.
measurable_outcome: Select an OpenClaw execution mode, provider path, and security boundary for a concrete local automation workflow within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# OpenClaw Agent Operations

## Workflow

1. Read `references/sources.md` and confirm the current OpenClaw repository and release state before implementation.
2. Decide whether the workload is truly local-first. If browser execution is the main requirement, compare against `Browser_Use_2026` first.
3. Define the security boundary explicitly: allowed directories, credentials, browser or session scope, and approval gates for shell or file actions.
4. Keep provider configuration and model selection separate from task routing.
5. Run one safe end-to-end task in a sandbox before broadening permissions.

## When to Use

- The system should run on the operator's machine.
- Browser work is only one part of a larger local automation surface.
- You need explicit control of local files, local tools, or local sessions.
- Privacy or operator control is more important than managed cloud browser infrastructure.

## Output Requirements

- State why OpenClaw is a better fit than Browser Use for the workload.
- State the primary provider path and one local security boundary.
- Include one sandbox or approval recommendation.
- Include one rollback or kill-switch note.
