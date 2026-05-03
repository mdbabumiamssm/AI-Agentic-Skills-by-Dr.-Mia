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

## Core Capabilities

- Reference governance control plane for AI coding agents: validate agent changes, enforce policy-as-code, track findings, and wire observability and compliance hooks. See controlkeel (Phoenix/Elixir) as a pattern when auditing agent-generated code changes in regulated environments.
- Live workflow-agent benchmarking: structure evals like Claw-Eval-Live by separating refreshable public workflow-demand signals (including ClawHub Top-500 skills in the current release) from reproducible, time-stamped release snapshots with fixed fixtures, services, workspaces, and graders. Grade with deterministic checks over execution traces, audit logs, service state, and post-run workspace artifact checks when evidence is sufficient, reserving structured LLM judging for semantic dimensions only. Treat leaderboard rank as insufficient: also report task-family breakdowns (e.g., HR, multi-system business workflows vs. local workspace repair) and overall completion alongside pass rate.
- Adversarial but legible task authoring for terminal/agent benchmarks: do not write tasks like prompts. Audit each task against a checklist for the common failure modes — AI-generated instructions, over-prescriptive specifications, clerical (non-conceptual) difficulty, oracle solutions that assume hidden knowledge, tests that validate the wrong thing, and reward-hackable environments. Difficulty must be conceptual rather than environmental, and verification logic should be reviewed adversarially before a task lands. Treat reward-hackability as the default failure: empirical review of popular terminal-agent benchmarks found >15% of tasks reward-hackable, so internal evals should run a hacking-pass review (intentional shortcut attempts, no-op solutions, exploiting graders) before they gate rollout.
- Eval task review checklist (use when authoring or auditing): (1) is the task adversarial — does it try to break the agent rather than help it succeed? (2) is difficulty conceptual, not clerical or environmental? (3) is the task legible — can a human confirm what success means without hidden context? (4) does the verifier check the actual goal, not a proxy artifact? (5) can the test be passed by reward-hacking, no-op exploits, or oracle leakage? (6) are instructions human-authored and minimal rather than over-prescriptive?

## References

- controlkeel — control plane for governed AI coding: https://github.com/aryaminus/controlkeel
- Claw-Eval-Live: A Live Agent Benchmark for Evolving Real-World Workflows: http://arxiv.org/abs/2604.28139v1
- What Makes a Good Terminal-Agent Benchmark Task — guideline for adversarial, difficult, and legible evaluation design: http://arxiv.org/abs/2604.28093v1
