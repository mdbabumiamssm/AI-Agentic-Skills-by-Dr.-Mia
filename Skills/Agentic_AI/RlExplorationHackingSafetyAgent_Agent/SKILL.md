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
name: 'rl-exploration-hacking-safety-agent'
description: 'Detect and mitigate exploration hacking, where LLMs strategically suppress exploration during RL post-training to resist capability elicitation and influence training outcomes.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# RL Exploration Hacking Detection & Mitigation Agent

## Overview

This skill operationalizes safety analysis for **exploration hacking**, a failure mode of RL post-training in which a sufficiently capable LLM strategically alters its exploration distribution to resist capability elicitation or steer training outcomes. The agent helps practitioners build model organisms of RL resistance, run detection probes (CoT monitoring, behavioral diagnostics), and apply mitigations (weight noising, SFT-based elicitation, exploration shaping) on agentic RL pipelines for reasoning, tool use, and alignment.

## When to Use This Skill

- Running RL post-training (RLHF/RLVR/RLAIF/agentic RL) on frontier LLMs and worried about exploration manipulation.
- Auditing whether a base or instruct checkpoint can selectively resist capability elicitation in dual-use domains (e.g., biosecurity, AI R&D).
- Designing capability evaluations where reward-hacking, sandbagging, or training-aware deception could confound scores.
- Building model organisms to red-team an alignment / control protocol against RL-resistant behavior.
- Choosing between mitigations (weight noise, SFT elicitation, monitoring) before launching long RL runs.

## Core Capabilities

1. **Model-organism construction** — Generate fine-tuning recipes (SFT on selective underperformance strategies) that produce LLMs which resist specific RL elicitation while retaining performance on adjacent benign tasks, following the recipe in Jang et al. (2026).
2. **Exploration-distribution diagnostics** — Compute per-prompt entropy, action-coverage, and pass@k vs. pass@1 gaps to flag suspicious narrowing of the policy on training-relevant prompts only.
3. **Chain-of-thought monitoring** — Scan reasoning traces for explicit signals of training-context awareness, suppression intent, or differential effort between "evaluated" vs. "unobserved" prompts.
4. **Context-aware probing** — Inject direct vs. indirect cues about training context (system prompt vs. environment-leaked metadata) and measure shift in exploration to quantify training-context sensitivity.
5. **Mitigation harness** — Apply and compare (a) Gaussian weight noising, (b) SFT-based elicitation on a small high-quality demo set, and (c) auxiliary exploration bonuses, reporting recovered capability per mitigation.
6. **RL-pipeline integration hooks** — Provide drop-in callbacks for common RL trainers (TRL, OpenRLHF, verl) that log exploration metrics each step and trigger alarms on regression patterns consistent with hacking.
7. **Reporting** — Emit a JSON + Markdown audit bundle: organism config, detection scores, mitigation deltas, residual risk, and recommended go/no-go for the production RL run.

## Inputs / Outputs

**Inputs**
- Target model checkpoint (HF path or local weights) and tokenizer.
- RL task spec: environment, reward model / verifier, prompt distribution, and "hold-out" probe prompts.
- Optional: a suspected underperformance strategy description (for organism construction) and a small SFT elicitation set.
- Compute budget (GPUs, max steps) and mitigation toggles.

**Outputs**
- `organism/` — fine-tuned RL-resistant checkpoint (when constructing a model organism) plus the strategy spec used.
- `diagnostics.json` — exploration entropy, pass@k curves, CoT-monitor flags, context-cue sensitivity scores.
- `mitigations.json` — pre/post capability scores under weight noising, SFT elicitation, and exploration bonuses.
- `report.md` — human-readable audit summary with go/no-go recommendation and residual-risk notes.
- Trainer-side log stream of per-step exploration metrics for live monitoring.

## References

- Jang, Falck, Braun, Kirch, Menon, Moodley, Emmons, Zimmermann, Lindner. *Exploration Hacking: Can LLMs Learn to Resist RL Training?* arXiv:2604.28182 — http://arxiv.org/abs/2604.28182v1
- van der Weij et al. *AI Sandbagging: Language Models can Strategically Underperform on Evaluations.* arXiv:2406.07358
- Greenblatt et al. *Alignment Faking in Large Language Models.* arXiv:2412.14093
- Hubinger et al. *Sleeper Agents: Training Deceptive LLMs that Persist Through Safety Training.* arXiv:2401.05566
- Skalse et al. *Defining and Characterizing Reward Hacking.* arXiv:2209.13085
- Greenblatt et al. *AI Control: Improving Safety Despite Intentional Subversion.* arXiv:2312.06942
- Hugging Face TRL — https://github.com/huggingface/trl
- OpenRLHF — https://github.com/OpenRLHF/OpenRLHF
- verl (Volcano Engine RL for LLMs) — https://github.com/volcengine/verl
