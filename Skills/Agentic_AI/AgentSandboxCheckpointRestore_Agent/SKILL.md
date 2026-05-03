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
name: 'agent-sandbox-checkpoint-restore'
description: 'Plan turn-aware, host-side checkpoint/restore for sandboxed agents by classifying OS-visible side effects and recovery-relevant state.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Agent Sandbox Checkpoint Restore

## Overview

Use this skill to design, review, or troubleshoot checkpoint/restore workflows for autonomous agents running inside containers, sandboxes, or microVMs. It focuses on the semantic gap identified by Crab: agent frameworks know turn and tool-call context, while the host operating system observes filesystem, process, and runtime side effects. The goal is to preserve recovery correctness without checkpointing every turn when no recovery-relevant state changed.

## When to Use This Skill

- Building fault tolerance for long-running agent jobs in sandboxed execution environments.
- Supporting rollback after risky tool calls, failed code edits, dependency changes, or destructive shell actions.
- Branching agent rollouts for reinforcement learning, evaluation, or alternate task trajectories.
- Running agents on spot, preemptible, or otherwise failure-prone infrastructure.
- Auditing whether chat-history-only recovery misses filesystem, process, or runtime side effects.
- Reducing checkpoint traffic for densely co-located agent sandboxes.

## Core Capabilities

1. **State surface mapping**: Identify which agent-visible actions can create OS-visible state across filesystems, processes, sockets, package managers, caches, and runtime artifacts.
2. **Turn-aware checkpoint policy**: Align checkpoint decisions with agent turn boundaries, tool calls, and LLM wait intervals instead of using fixed per-turn or purely time-based snapshots.
3. **Recovery relevance classification**: Separate no-op or transient turns from turns that alter state needed for correct restart, rollback, or rollout branching.
4. **Host-side integration plan**: Prefer transparent host-level observation and coordination when the agent framework or workload cannot be modified.
5. **Checkpoint scheduling**: Account for co-located sandboxes by limiting bursty checkpoint traffic and overlapping checkpoint work with agent idle or LLM wait time.
6. **Restore validation**: Define replay tests that compare restored execution against expected filesystem state, process/runtime continuity, and agent task progress.

## Inputs / Outputs

**Inputs**

- Agent runtime model: container, sandbox, microVM, or host process layout.
- Tool-call and turn structure, including where LLM waits and shell/tool boundaries occur.
- Recovery objective: fault tolerance, rollback, rollout branching, spot execution, or forensic replay.
- State inventory: files, processes, package installs, temporary directories, sockets, caches, environment variables, and external mounts.
- Existing checkpoint/restore backend constraints, such as snapshot granularity, pause time, storage cost, and restore semantics.

**Outputs**

- A turn-aware checkpoint/restore design or review memo.
- A classification table for recovery-relevant and non-recovery-relevant side effects.
- A checkpoint policy specifying when to skip, take lightweight metadata, or capture full state.
- A restore test plan covering correctness, rollback boundaries, and failure injection.
- Operational guidance for scheduling checkpoint traffic across multiple active sandboxes.

## References

- Crab: A Semantics-Aware Checkpoint/Restore Runtime for Agent Sandboxes: http://arxiv.org/abs/2604.28138v1
