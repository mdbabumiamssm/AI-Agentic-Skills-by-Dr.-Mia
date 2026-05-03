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
name: 'synthetic-computers-long-horizon-training'
description: 'Generate persona-conditioned synthetic computer environments and run multi-thousand-turn long-horizon productivity simulations to produce experiential training data for agent self-improvement.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Synthetic Computers Long-Horizon Productivity Agent Training

## Overview

This skill operationalizes the *Synthetic Computers at Scale* methodology: it constructs realistic, persona-conditioned synthetic computer environments (folder hierarchies, documents, spreadsheets, presentations) and then runs long-horizon (>2,000-turn, ~8+ hour) productivity simulations against them. The resulting trajectories form a substrate of experiential learning signals usable for agent supervised fine-tuning and agentic reinforcement learning on month-scale professional work scenarios.

## When to Use This Skill

- Building large-scale training corpora for productivity agents that must operate over realistic filesystems.
- Generating month-of-work-equivalent rollouts requiring multi-deliverable coordination.
- Stress-testing an agent on long-horizon planning, filesystem grounding, and collaborator coordination.
- Producing in-domain and out-of-domain productivity evaluation environments.
- Scaling synthetic-environment data generation across diverse personas, professions, and roles.
- Creating self-improvement loops for agentic RL where the environment itself is procedurally synthesized.

## Core Capabilities

1. **Persona-Conditioned Computer Synthesis** — Sample a persona (profession, role, context) and instantiate a synthetic computer with a coherent directory tree and content-rich artifacts (documents, sheets, slide decks) consistent with that user's working life.
2. **Objective Generation Agent** — Drives a planner agent that authors month-scale productivity objectives tied to the synthesized computer, requiring multiple professional deliverables grounded in existing artifacts.
3. **Long-Horizon User Simulation** — Runs an actor agent that plays the synthetic user, navigating the filesystem, editing artifacts, and coordinating with simulated collaborators across thousands of turns until objectives are completed.
4. **Collaborator Modeling** — Spawns auxiliary agents representing teammates, managers, or external contacts to inject realistic communication and review loops.
5. **Trajectory Capture for Training** — Logs full action/observation traces, intermediate artifacts, and outcomes in a format suitable for SFT, preference learning, or agentic RL.
6. **Scaling Harness** — Parallelizes synthesis and rollout across many synthetic computers (the source paper reports 1,000 computers, with billion-persona scaling in principle).
7. **Evaluation Hooks** — Holds out synthetic computers as in-domain and out-of-domain productivity benchmarks to quantify agent improvement after training on collected rollouts.

## Inputs / Outputs

**Inputs**
- A persona specification (profession, role, seniority, context, productivity needs).
- Optional artifact templates and a content-generation backbone model.
- Compute budget (rollout horizon cap, parallel-environment count, max wall-clock).
- A base agent policy to drive both objective generation and user simulation.

**Outputs**
- A populated synthetic computer: directory tree plus document/spreadsheet/presentation artifacts.
- A set of long-horizon objectives bound to that computer.
- Full multi-thousand-turn rollout trajectories with actions, observations, and produced deliverables.
- A training-ready dataset of trajectories for SFT or agentic RL.
- Evaluation splits for in-domain and out-of-domain productivity scoring.

## References

- Ge, T.; Peng, B.; Cheng, H.; Gao, J. *Synthetic Computers at Scale for Long-Horizon Productivity Simulation.* arXiv:2604.28181v1 — http://arxiv.org/abs/2604.28181v1
