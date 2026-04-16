---
name: agent-evals-observability-2026
description: Evaluate and observe agent systems with traced runs, curated datasets, rollout gates, and regression loops. Use when an agent workflow is moving beyond prototype and needs measurable reliability before broader autonomy or production rollout.
measurable_outcome: Define an eval set, trace sink, rollout gate, and regression loop for a concrete agent workflow within 30 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Agent Evals and Observability (2026)

## Workflow

1. Read `references/sources.md` and choose the primary tracing and evaluation surface before changing prompts or tools.
2. Define the evaluation unit explicitly: final answer quality, tool selection, tool arguments, trajectory quality, latency, cost, or safety.
3. Build a compact benchmark set first from golden tasks, known failures, and a small live-traffic sample.
4. Add traces and metadata so failed runs can be reproduced with model, prompt, tool, and environment context.
5. Set rollout gates before enabling side effects: minimum pass rate, maximum regression, latency budget, and manual review triggers.
6. Keep the loop closed by promoting failing production traces back into the offline dataset.

## When to Use

- You are preparing an agent workflow for production or wider internal use.
- Model or framework changes need regression coverage.
- Tool calling, approvals, or browser execution introduce failure modes that prompt-only tests will miss.
- You need both pre-deployment evals and post-deployment monitoring.

## Output Requirements

- Return the evaluation unit and why it matters for the workload.
- State the primary trace sink and the minimum metadata that must be captured.
- Include one offline eval path and one online monitoring path.
- Include one rollout gate and one rollback trigger.
