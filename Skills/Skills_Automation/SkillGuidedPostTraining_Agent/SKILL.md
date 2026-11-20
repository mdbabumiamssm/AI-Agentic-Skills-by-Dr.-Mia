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
name: 'skill-guided-post-training'
description: 'Design and evaluate natural-language skill interventions that improve teacher demonstrations for supervised fine-tuning, measuring task success separately from undesirable behavioral habits.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Skill-Guided Post-Training

## Overview

Use natural-language skills to shape teacher demonstrations before supervised fine-tuning (SFT), because successful trajectories can still demonstrate undesirable habits. AutoMaster develops interventions through rollout review and student feedback; Skill-FT retains generating skill context during training and evaluates students without generated skills ([paper, Sections 3–5](https://arxiv.org/html/2611.26510v1)). This skill supplies a practical research protocol; its audit rubrics and acceptance gates below are local recommendations, not reported experimental procedures.

## When to Use This Skill

- Successful teacher traces contain unsupported assumptions, repeated ineffective actions, or weak verification.
- Useful behaviors appear too rarely in demonstrations selected for SFT.
- You need to compare guided demonstrations against successful unguided demonstrations under controlled conditions.
- You want to test whether retaining skill context during SFT improves behavior without supplying skills at inference.

## Core Capabilities

1. **Verify evidence and execution readiness.** Read the versioned paper, especially Section 3 and Appendix A, before choosing implementation details. Record source version, access date, artifact revisions, licenses, and missing dependencies. As checked on 2026-11-20, the [official repository](https://github.com/RUCAIBox/AutoMaster) describes nine selected skills, inference recipes, and an AweAgent submodule; the [model collection](https://huggingface.co/collections/RUC-AIBOX/automaster) lists six checkpoints. These listings do not establish availability of the complete autonomous training controller or training corpus. Inspect actual artifacts before prescribing commands, interfaces, or hardware; mark unavailable components explicitly and produce an adapted protocol instead of claiming reproduction.

2. **Review outcomes and behavior independently.** Sample both successful and failed traces across task types, including initial student traces where available. Use the task verifier for outcome labels and an independently defined rubric for behavior. Record an evidence span or tool event for each judgment; use `unknown` when the trace cannot establish a claim. Do not treat a confident final answer or a self-reported check as evidence of correctness.

   | Dimension | Operational audit question | Suggested measurement |
   | --- | --- | --- |
   | Task success | Did the external verifier accept the result? | Pass/fail or native task score |
   | Unsupported assumptions | Was a consequential premise used without evidence or an explicit test? | Affected traces / reviewed traces |
   | Unproductive repetition | Did the agent repeat an ineffective action without new information? | Events per trace, with context |
   | Superficial verification | Did the check exercise the actual requirement and inspect its result? | Inadequate checks / applicable checks |
   | Effective behavior | Did a relevant opportunity trigger the intended action? | Demonstrated behaviors / eligible opportunities |

   Define denominators before comparison. Review ambiguous cases consistently and record reviewer disagreement rather than forcing agreement. Report the subset where both guided and unguided teachers succeed separately, so success improvements cannot conceal persistent habits.

3. **Design scoped interventions.** Convert each observed weakness into a falsifiable hypothesis. Specify the trigger, preferred action, evidence needed, stopping condition, and cases where the instruction should not apply. Assign a version and targeted task subset; avoid task answers or evaluation-specific clues in reusable skills. Example local intervention: after a repeated search yields no new evidence, identify the unsupported premise, choose a different discriminating clue, and verify the resulting candidate. Predict both a behavioral improvement and a possible regression such as excessive checking.

4. **Generate and compare teacher rollouts.** Pilot guided and unguided runs on matched task IDs with fixed teacher, environment, decoding settings, and resource limits. Preserve complete messages, tool results, verifier outputs, skill versions, and failure reasons. The paper appends skills to the base system prompt and scales promising pilots before selecting demonstrations ([Section 3](https://arxiv.org/html/2611.26510v1)). For this workflow, predefine acceptable outcome regression and behavioral improvement thresholds; revise or reject interventions that fail them. Log timeouts and invalid outputs instead of silently excluding them.

5. **Select traceable demonstrations.** Write an inclusion decision and reason for each candidate. Prefer evidence of useful behavior and valid task execution; keep any exceptions visible. Deduplicate and check train/evaluation overlap before training. Preserve useful unguided traces where warranted, but record the guided/unguided mixture. Fix training-set size across comparison arms and report token totals, task composition, and exclusion rates to expose remaining confounds. If too few eligible demonstrations exist, report the shortfall rather than silently duplicating examples or changing the comparison budget.

6. **Run controlled SFT comparisons.** Freeze a manifest containing the initial student revision, teacher revision, task pool, selected trace IDs, tokenizer/chat template, training configuration, and evaluation harness. The paper restarts each outer round from the same initial student checkpoint; assistant outputs receive loss while system, user, and tool messages provide context ([Section 3 and Appendix A.3](https://arxiv.org/html/2611.26510v1)). Validate message serialization, loss masking, and truncation on a small batch before launching a supported training job. Use the verified implementation's configuration rather than invented APIs. Record any deviation from the paper as an adaptation.

   | Comparison arm | Purpose |
   | --- | --- |
   | Untrained student | Establish the starting point |
   | Vanilla SFT on successful unguided traces | Main baseline with the same example count and recipe |
   | Guided and selected traces, skill context retained | Main intervention |
   | Identical guided traces, skill context removed | Isolate prompt context; keep supervised outputs identical |
   | Optional selection-only and guidance-only arms | Separate generation and selection effects |
   | Optional task-matched unguided arm | Match task IDs and per-task counts; disclose outcome differences |

7. **Evaluate transfer and iterate within budget.** Evaluate every student without generated skills using identical settings. Report task scores and the behavior rubric separately, plus cost and completion failures. Use repeated evaluations when feasible; distinguish variability of a fixed checkpoint from variability across independent training runs. Compare paired task results and state uncertainty rather than declaring a gain from one favorable run. As a local safeguard against adaptive evaluation overfitting, use development feedback for revisions and reserve an untouched final test set. Link every result to its checkpoint, manifest, and intervention version. Stop at the declared compute/time budget or acceptance criterion; record regressions and inconclusive findings as valid research outcomes.

## Inputs / Outputs

**Inputs:** A task pool with split identifiers and verifiers; teacher and initial student identifiers; available trajectories with tool observations; candidate skills or behavior hypotheses; a supported rollout/training/evaluation stack; and explicit time, compute, and data-use constraints. If an input is missing, identify which stage it blocks while continuing independent review.

**Outputs:** An evidence and artifact inventory; an annotated trajectory audit; versioned interventions with applicability rules; a rollout and selection manifest; a controlled SFT/evaluation plan; and, when execution is available, linked checkpoints and a report separating outcome changes, behavior changes, uncertainty, and cost. Label every result as observed, planned, or unavailable.

**15-minute completion contract:** Complete a bounded research setup pass, not an entire model-training campaign. Use the first portion to verify sources and assets, then audit a small declared trace sample, design one intervention, and specify comparison arms and acceptance gates. A valid output contains provenance, separate outcome/behavior criteria, at least one testable intervention, fixed experimental controls, and either pilot evidence or the precise blocker and next executable step. If no traces are available, deliver an explicitly untested protocol; never fabricate an audit or claim completed training. Longer rollout and training jobs require their own stated budgets and completion status.

## References

- [AutoMaster: Autonomous Post-Training by Making the Same Model a Better Teacher, arXiv:2611.26510v1](http://arxiv.org/abs/2611.26510v1) — source finding; [full method and appendices](https://arxiv.org/html/2611.26510v1).
- [Official AutoMaster repository](https://github.com/RUCAIBox/AutoMaster) — released skill examples, inference recipes, and harness linkage; recheck release scope before use.
- [AutoMaster model collection](https://huggingface.co/collections/RUC-AIBOX/automaster) — checkpoint discovery; inspect individual model cards and revisions before execution.
