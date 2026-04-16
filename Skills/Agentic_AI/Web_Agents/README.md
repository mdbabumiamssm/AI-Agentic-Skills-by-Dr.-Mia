# Web Agents (2026)

This directory covers browser-native and browser-adjacent agent execution.

## Use Browser Use first when:

- live browser execution is the core requirement
- managed browser infrastructure or cloud workspaces matter
- browser sessions need to be portable under a larger orchestrator

## Use OpenClaw first when:

- the workload is local-first
- browser work is part of a broader local automation surface
- local files, shells, and operator-controlled sessions are central

## Selection rule

Do not treat every web-facing task as a browser-agent problem.

- If the real need is orchestration, choose an orchestration framework first.
- If the real need is local automation, use OpenClaw.
- If the real need is managed browser execution, use Browser Use.
