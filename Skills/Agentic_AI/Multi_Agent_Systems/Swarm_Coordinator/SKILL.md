---
name: swarm-coordinator-pattern
description: Design a supervisor-worker or debate-style multi-agent coordinator with explicit delegation and synthesis boundaries. Use when a task genuinely benefits from controlled specialization rather than a single-agent loop.
measurable_outcome: Define a coordinator pattern, delegation boundary, and synthesis path for a concrete multi-agent workflow within 20 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Swarm Coordinator Pattern

## Workflow

1. Start from the parent `Multi_Agent_Systems/` guidance before adding more agents.
2. Prove that the workload needs multiple specialists rather than a single agent with tools.
3. Define agent roles, handoff boundaries, and the final synthesis contract explicitly.
4. Add guardrails for duplicate work, deadlocks, and conflicting outputs.
5. Keep the coordinator observable: each delegation should be traceable and reviewable.

## Output Requirements

- Return the chosen multi-agent pattern and why it is justified.
- State the coordinator responsibilities and worker boundaries.
- Include one conflict-resolution rule.
- Include one observability or evaluation note.
