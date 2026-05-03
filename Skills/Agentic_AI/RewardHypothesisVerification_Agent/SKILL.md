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
name: 'reward-hypothesis-verification'
description: 'Verify and deploy LLM-generated reward hypotheses using competence-aware short-horizon forks and phase-aware RL training decisions.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Reward Hypothesis Verification

## Overview

Use this skill to treat LLM-generated reward functions as reward hypotheses that must be verified before they are trusted as reinforcement learning objectives. It guides competence-aware short-horizon fork verification, phase-aware deployment, schedule selection, and boundary testing so reward candidates are evaluated under the policy competence where they will actually be used.

## When to Use This Skill

- An RL workflow uses LLM-generated, evolved, or selected reward candidates.
- Reward rankings appear unstable across training checkpoints or policy competence levels.
- A reward candidate performs well in one training phase but fails, stalls, or degrades later.
- A team needs a deployment protocol for switching, delaying, or rejecting generated rewards.
- Sparse-reward, manipulation, robotics, or agent-training tasks need verification before reward deployment.
- Reward evaluation must include conservative baselines, held-out schedule selection, compute-matched controls, or failure-boundary tests.

## Core Capabilities

1. **Frame generated rewards as hypotheses** - Treat each candidate reward as a provisional training objective whose utility can depend on policy competence, environment phase, and checkpoint state.
2. **Run shared-checkpoint fork verification** - Compare small sets of reward hypotheses by forking from the same policy checkpoint and training or rolling out for a short horizon under matched compute.
3. **Estimate competence thresholds** - Identify policy competence ranges where reward rankings become informative enough to guide deployment, while marking low-competence rankings as unreliable unless supported by task evidence.
4. **Select phase-aware deployment schedules** - Choose when to introduce, switch, retain, or reject reward hypotheses using held-out verification outcomes rather than assuming a fixed warm-up schedule.
5. **Use conservative selector baselines** - Compare proposed reward deployment against simple selectors, fixed schedules, original rewards, and compute-matched controls to separate verification value from extra optimization budget.
6. **Test failure boundaries** - Include dense-reward, all-failure, and low-signal cases to determine when verification cannot distinguish candidates or when generated rewards should not be deployed.
7. **Report deployment evidence** - Produce a compact decision record containing checkpoints tested, competence estimates, fork outcomes, selected schedule, rejected candidates, controls, and residual risks.

## Inputs / Outputs

**Inputs**

- Task or environment description, including success criteria and available observations.
- Candidate reward hypotheses, reward code, or natural-language reward specifications.
- Policy checkpoints, training logs, competence estimates, and evaluation metrics.
- Verification budget, fork horizon, random seeds, and compute constraints.
- Baselines such as original reward, fixed schedule, conservative selector, or no-switch control.

**Outputs**

- Ranked reward hypotheses for each tested competence phase, with uncertainty or instability notes.
- Recommended deployment action: deploy now, delay until threshold, switch at a phase boundary, retain current reward, or reject candidates.
- Held-out schedule selection rationale and compute-matched comparison summary.
- Boundary-test findings describing cases where the protocol is uninformative or unsafe to rely on.
- Reproducible verification record with checkpoint IDs, fork settings, metrics, and final deployment recommendation.

## References

- RHyVE: Competence-Aware Verification and Phase-Aware Deployment for LLM-Generated Reward Hypotheses - http://arxiv.org/abs/2604.28056v1
- Eureka: Human-Level Reward Design via Coding Large Language Models - https://arxiv.org/abs/2310.12931
- Reinforcement Learning from Human Feedback: Learning Dynamic Choices via Pessimism - https://arxiv.org/abs/2305.18438
- RewardBench: Evaluating Reward Models for Language Modeling - https://arxiv.org/abs/2403.13787
