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
name: 'skill-guided-agent-post-training'
description: 'Design and assess skill-guided teacher trajectories for agent fine-tuning, with controlled training-data comparisons and held-out evaluation of task outcomes and behavior.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Skill-Guided Agent Post-Training

## Overview

Turn observed agent habits into testable natural-language interventions and a controlled post-training research protocol. Skill-FT shapes teacher demonstrations, while AutoMaster iterates on skills and data selection using rollout and student feedback ([method, Section 3](https://arxiv.org/html/2611.26510v1)). The audit rubric, split discipline, and completion criteria below are proposed research practices, not claims about the paper's implementation.

## When to Use This Skill

- Successful trajectories still contain unsupported assumptions, ineffective repetition, or weak verification.
- You need operational specifications for behaviors that are rare in existing demonstrations.
- You want to compare skill-guided supervision with successful unguided demonstrations.
- You need to assess behavioral transfer without supplying generated skills at evaluation.
- You are preparing a bounded pilot before committing to a larger training campaign.

## Core Capabilities

1. **Establish evidence and scope.** Read the versioned method and Appendix A before adopting implementation details. Record accessible artifacts, model revisions, harness versions, budget, and missing dependencies; distinguish reproduction from adaptation. Check the populated skill inventory for overlap. A related entry exists at `Skills/Skills_Automation/SkillGuidedPostTraining_Agent/SKILL.md`; coordinate future maintenance without modifying it as a side effect of this workflow.

2. **Audit behavior independently of success.** Review a declared sample containing successes and failures. Attach trace IDs and observable message or tool-event evidence to every judgment; use `unknown` for missing evidence. Apply this proposed rubric consistently across comparison arms:

   | Dimension | Evidence to inspect | Measurement |
   | --- | --- | --- |
   | Task outcome | External verifier or native task scoring | Score and completion status |
   | Unsupported assumption | Consequential premise used without support or a discriminating check | Affected traces / assessable traces |
   | Repetition | Repeated ineffective action without new evidence or a changed hypothesis | Events per trace |
   | Verification quality | Check addresses the requirement and its result is inspected | Adequate checks / applicable opportunities |
   | Desired behavior | Trigger occurs and the specified action follows | Satisfied triggers / eligible triggers |

   Report denominators and ambiguous cases. Inspect pairs where both conditions succeed separately, so improved outcomes cannot hide persistent habits.

3. **Specify falsifiable interventions.** For each candidate, record `skill_id`, version, evidence motivating it, trigger, action, required evidence, stopping condition, scope exclusions, and predicted regression. Keep answers and evaluation-specific clues out of reusable instructions. For example, propose that a stalled investigation should identify an untested premise and choose a check that distinguishes competing explanations; test whether this reduces redundant actions without excessive tool use.

4. **Pilot and select demonstrations.** The method appends skills to the teacher's base system prompt, pilots interventions, then selects useful guided or unguided trajectories. For this workflow, compare matched tasks with fixed decoding and execution budgets, preserve unsuccessful attempts, and predeclare acceptable regressions. Record each inclusion or exclusion with its reason; deduplicate traces and audit split overlap. Report task composition and token totals as potential confounds even when example counts match. A shortage of eligible data is a reported limitation, not a reason to silently duplicate examples.

5. **Construct controlled training inputs.** The paper retains generating skill context, supervises assistant outputs while masking system/user/tool messages, and restarts each round from the original student checkpoint ([Section 3 and Appendix A.3](https://arxiv.org/html/2611.26510v1)). Validate serialization, role boundaries, loss masks, and truncation before execution. Preserve tool observations as context; never invent missing trajectory content. Use verified training interfaces and record departures from the published recipe.

   | Arm | Control or question |
   | --- | --- |
   | Initial student | Starting performance before additional training |
   | Successful unguided teacher demonstrations | Vanilla SFT baseline |
   | Selected skill-guided demonstrations with skill context | Main intervention |
   | Same selected demonstrations with only skill text removed | Context ablation with identical supervised outputs |

   Hold teacher, task pool, training-set size, initial checkpoint, and SFT recipe fixed across trained arms. Optionally add selection-only or task-matched comparisons, documenting their selection rules and any differences in outcome composition. Record trace IDs, prompt versions, tokenizer, training configuration, and checkpoint IDs in a manifest.

6. **Evaluate transfer and revise.** The method evaluates students without generated skills and uses scores plus trajectories for subsequent revisions. As a local research safeguard, separate development feedback from an untouched final test set; split related repositories or task variants together. Report outcome scores, behavioral rates, resource use, invalid runs, and uncertainty separately. Distinguish repeated evaluation of one checkpoint from independent training runs. Diagnose whether a weakness reflects ineffective guidance, excluded demonstrations, or failed transfer before changing the intervention. Stop at the declared budget or acceptance criterion, retaining negative and inconclusive results.

## Inputs / Outputs

**Inputs:** Task pool and split identifiers; verifier definitions; teacher and initial student revisions; trajectories with tool observations and outcomes; behavior hypotheses or candidate skills; supported rollout, training, and evaluation interfaces; and time/compute constraints.

**Outputs:** A source and dependency record, evidence-linked trajectory audit, versioned intervention specification, data-selection manifest, controlled experiment plan, and an evaluation report when runs exist. Link reported measurements to actual run artifacts and label each stage `planned`, `completed`, or `blocked`; checkpoints are outputs only after training has completed.

**15-minute completion contract:** Deliver a bounded research setup or pilot review, not a promise to finish model training. Declare the reviewed sample size, specify at least one intervention and its measurable acceptance criteria, define comparison arms and split checks, and identify the next executable step. With no usable trajectories or execution stack, deliver an explicitly untested protocol and list the missing inputs; do not fabricate audit results or claim behavioral improvement. Longer training and evaluation require separate budgets.

**Validity check:** Ensure every audit judgment has evidence or an explicit unknown label; every selected trace has provenance; context-ablation outputs remain identical; final evaluation is isolated from intervention design; and reported gains come from completed comparisons rather than expectations.

## References

- [AutoMaster: Autonomous Post-Training by Making the Same Model a Better Teacher, arXiv:2611.26510v1](http://arxiv.org/abs/2611.26510v1) — source finding.
- [Full method and experimental details](https://arxiv.org/html/2611.26510v1) — Section 3 for intervention loops; Appendix A for rollout, SFT, and evaluation procedures.
