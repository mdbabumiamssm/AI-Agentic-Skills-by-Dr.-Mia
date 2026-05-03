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
name: 'latent-adversarial-multiturn-detector'
description: 'Detect covert multi-turn prompt injection by probing residual-stream activations for "adversarial restlessness" trajectory signatures that text-level filters miss.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Latent Activation Multi-Turn Attack Detector

## Overview

This skill operationalizes activation-level detection of multi-turn prompt injection attacks against deployed LLM agents. It captures the residual-stream trajectory of a conversation, computes scalar "adversarial restlessness" features (path length, phase-shift magnitude, directional drift), and trains/serves a lightweight per-model probe that flags trust-building → pivoting → escalation attack paths even when individual turns appear benign. The signal is complementary to text-level guardrails and aimed at production agent deployments where covert manipulation is in scope.

## When to Use This Skill

- You are hardening a multi-turn agent (chatbot, copilot, tool-using agent) against jailbreaks that span several turns.
- Text-only safety classifiers, regex filters, or per-turn moderation pass adversarial dialogues that escalate gradually.
- You have white-box or grey-box access to the target model's residual-stream activations (open-weights models 24B–70B, or hosted models exposing hidden states).
- You need a conversation-level decision (escalate / block / hand-off) rather than per-turn moderation.
- You are evaluating defense generalization across attack distributions (synthetic, real-world chat logs, dialog safety benchmarks).

## Core Capabilities

1. **Activation capture pipeline** — Hooks the residual stream at a configurable layer of a HuggingFace transformer, recording per-turn activation vectors for an entire conversation with stable offsets and turn boundaries.
2. **Trajectory feature extraction** — Computes the five scalar features that characterize adversarial restlessness: total path length, mean inter-turn step size, max step size, directional variance, and cumulative drift from the conversation centroid.
3. **Per-model probe training** — Fits a lightweight classifier (logistic regression / shallow MLP) on labeled conversations using three-phase turn-level labels (benign / pivoting / adversarial) which the source paper shows are necessary to avoid 50–59% false-positive rates from binary labels alone.
4. **Multi-source training harness** — Combines synthetic, LMSYS-Chat-1M, and SafeDialBench-style sources with leave-one-source-out evaluation so deployers can audit which attack distributions are actually represented.
5. **Inference-time guardrail** — Streams trajectory features as the conversation grows, emitting a calibrated probability and a configurable action (warn, require human review, terminate) at a chosen false-positive operating point.
6. **Architecture-aware caveats** — Probes are explicitly model-specific; the skill enforces that a probe trained on one base model is not silently reused on another family.

## Inputs / Outputs

**Inputs**
- Target model handle (HF model id or local checkpoint) with hidden-state access.
- Layer index (or list of indices) to probe.
- Labeled conversation dataset with three-phase turn labels (benign / pivoting / adversarial) for training; raw conversation transcript for inference.
- Operating-point configuration (target false-positive rate, decision threshold).

**Outputs**
- Cached activation tensors and per-conversation trajectory feature CSVs.
- Trained probe artifact (weights + scaler + metadata: base model, layer, sources, FPR/TPR table).
- Inference verdict per conversation: `{score, label, contributing_features, turn_at_which_signal_crossed_threshold}`.
- Evaluation report with per-source detection rates and leave-one-source-out generalization numbers.

## References

- Source paper: Kulkarni, P. *Latent Adversarial Detection: Adaptive Probing of LLM Activations for Multi-Turn Attack Detection.* arXiv:2604.28129 — http://arxiv.org/abs/2604.28129v1
- LMSYS-Chat-1M dataset: https://huggingface.co/datasets/lmsys/lmsys-chat-1m
- SafeDialBench: https://arxiv.org/abs/2502.11090
- Background on linear probes over residual streams: Alain & Bengio, *Understanding intermediate layers using linear classifier probes*, https://arxiv.org/abs/1610.01644
- Multi-turn jailbreak threat model (Crescendo attack): https://arxiv.org/abs/2404.01833
