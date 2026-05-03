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
name: 'llm-rl-posttraining-optimization'
description: 'Optimize LLM RL post-training for efficient reasoning with GRPO variants, latent constraints, advantage estimation, cost controls, and verification gates.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# LLM RL Posttraining Optimization

## Overview

Use this skill to design, audit, or troubleshoot reinforcement-learning post-training workflows for LLM reasoning systems. It emphasizes stable GRPO-style optimization, latent reasoning constraints, low-sample advantage estimation, cost-aware rollout design, and verification gates that make results efficient and auditable.

## When to Use This Skill

- A task asks for RL post-training of reasoning, math, code, tool-use, or multimodal models.
- A workflow uses GRPO, PPO, DPO-adjacent evaluation, reward-model-free RL, or group-relative advantage updates.
- Latent reasoning, hidden-state reasoning, compressed chain-of-thought, or shortened reasoning traces are part of the design.
- Training is unstable because rewards are sparse, rollouts leave the valid solution space, or correct trajectories conflict.
- The user needs cost controls for rollout count, pass@k evaluation, sampling temperature, verifier calls, or reward computation.
- The user needs an audit checklist for reward leakage, invalid samples, off-manifold latent states, or optimization-verification mismatch.

## Core Capabilities

1. **GRPO Variant Selection** - Decide whether standard GRPO, latent-aware GRPO, cost-aware GRPO, or a conservative PPO-style baseline best fits the task, model, and reward signal.
2. **Latent Reasoning Stabilization** - Check for invalid latent samples, off-manifold exploration, mixture non-closure, and reward updates that reinforce averaged or invalid hidden states.
3. **Advantage Estimation Design** - Choose group-relative, masked, kernelized, or low-sample advantage estimators based on rollout diversity, reward sparsity, and variance risk.
4. **Exploration-Optimization Alignment** - Ensure token-level or latent-step updates correspond to the actual successful path rather than indiscriminately crediting all sampled steps.
5. **Cost-Aware Rollout Planning** - Set rollout budgets, verifier budgets, sampling strategies, and early-stop criteria so training remains measurable and reproducible.
6. **Verification Gates** - Add checks for reward correctness, pass@1/pass@k reporting, invalid sample masking, prompt contamination, benchmark leakage, and regression against explicit-reasoning baselines.
7. **Operational Audit Trail** - Produce concise artifacts describing objectives, reward definitions, sampling policy, acceptance criteria, failure modes, and rollback triggers.

## Inputs / Outputs

**Inputs**

- Model family, checkpoint, tokenizer, and whether reasoning is explicit, latent, or hybrid.
- Task distribution, benchmark names, reward function, verifier design, and success criteria.
- Rollout configuration: group size, sampling method, temperature, max tokens or latent steps, and pass@k settings.
- Optimization configuration: objective, advantage estimator, masking rules, KL or trust-region controls, batch size, and compute budget.
- Existing logs or symptoms: reward curves, invalid sample rates, verifier disagreement, mode collapse, or degraded reasoning length/quality tradeoffs.

**Outputs**

- A recommended RL post-training recipe with objective, rollout plan, advantage estimator, masking strategy, and verification gates.
- A risk checklist covering latent manifold validity, reward leakage, optimization misalignment, cost overruns, and benchmark contamination.
- A minimal experiment matrix comparing latent and explicit reasoning baselines under the same evaluation protocol.
- A reproducibility note listing seeds, sampling settings, model checkpoints, verifier versions, and acceptance thresholds.

## References

- Latent-GRPO: Group Relative Policy Optimization for Latent Reasoning - http://arxiv.org/abs/2604.27998v1
- DeepSeekMath: Pushing the Limits of Mathematical Reasoning in Open Language Models - https://arxiv.org/abs/2402.03300
- DeepSeek-R1: Incentivizing Reasoning Capability in LLMs via Reinforcement Learning - https://arxiv.org/abs/2501.12948
- Proximal Policy Optimization Algorithms - https://arxiv.org/abs/1707.06347
- Training language models to follow instructions with human feedback - https://arxiv.org/abs/2203.02155
