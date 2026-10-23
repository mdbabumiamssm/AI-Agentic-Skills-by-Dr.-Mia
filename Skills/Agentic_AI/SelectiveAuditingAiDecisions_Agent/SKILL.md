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
name: 'selective-auditing-ai-decisions'
description: 'Monitor sequential AI decisions using confidence-guided expert audits, record corrections, measure audit coverage, and document the limits of error-control guarantees.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Selective Auditing of AI Decisions

## Overview

Use this skill to organize online review of AI proposals accompanied by confidence scores, including scores from a separate verifier. It produces an auditable monitoring policy and evidence report aimed at reducing expert workload while tracking errors that escape review.

## When to Use This Skill

- Route a stream of AI decisions between automatic acceptance and expert review.
- Monitor proposals whose confidence comes from a separate verifier.
- Evaluate selective auditing for safety screening, qualitative annotation, or labeling.
- Diagnose missing labels, biased audit samples, or unsupported reliability claims.
- Document whether a monitoring implementation satisfies a cited algorithm's conditions.

## Core Capabilities

1. **Define the monitoring contract.** Specify the unit of decision, binary error rubric, user-selected tolerance `alpha`, review timing, expert availability, and actions permitted before review. Distinguish an incorrect original proposal from its corrected final decision. Do not infer an acceptable error tolerance from an audit budget.

2. **Record decisions before routing.** Maintain one immutable proposal record per round with `decision_id`, timestamp, proposal or access-controlled reference, proposer version, verifier version, confidence, score direction, policy version, pre-update state, audit probability, random draw or replay identifier, and routing result. Append review status, original-proposal error label, correction, final decision, and review completion time. Keep unavailable labels null; unaudited does not mean correct.

3. **Select a reproducible policy.** The following compact specification summarizes Section 3.1 and Theorem 1 of the [verified full text](https://arxiv.org/html/2610.31662v1):

   - Let `C` be confidence in `[0,1]`, `e` indicate proposal error, and `H` indicate auditing. Initialize `tau` in `[0,1]`; choose `eta > 0`, `0 < alpha < 1`, and `0 < gamma <= 1-alpha`.
   - Set `pi = 1` if `C < tau`, otherwise `gamma`; draw `H ~ Bernoulli(pi)` using fresh conditional randomness.
   - If audited, obtain `e`, correct errors, and update `tau += eta*(1-pi-alpha)*e/pi`. Otherwise leave `tau` unchanged. Do not clip `tau`.
   - For `n = E[sum(e)] > 0`, the bound is `|E[sum(e*(1-H))]/n - alpha| <= (1/eta + 1/gamma)/n`.
   - This is a ratio of expectations, not a per-run or high-probability cap. It assumes correct audits and immediate feedback; sparse errors weaken the bound. Basic thresholding is not optimally efficient.

   Record parameter choices rather than inventing recommended values. If using calibrated scores, batches, or the optimal method, retrieve the corresponding algorithm and theorem before implementation; do not transplant the basic update's guarantee to a modified policy.

4. **Handle review failures explicitly.** A selected but unfinished audit is pending, not completed coverage. Hold its decision until review if the monitoring contract requires correction before release. For missing scores, unavailable experts, disputed labels, or delayed feedback, record the exception and use the agreed hold or review fallback. Document where the implemented process departs from the verified method; do not silently label such intervals as guaranteed.

5. **Measure coverage without hiding missing outcomes.** Report proposed, selected, completed, pending, corrected, and automatically released counts. Calculate selection rate as selected audits divided by proposals and completed coverage as completed audits divided by proposals; report latency and realized review cost separately. Report audited error prevalence with its selected-sample denominator, not as population error prevalence.

   With complete retrospective labels, calculate escaped-error fraction as erroneous unaudited proposals divided by all erroneous proposals, and residual error per decision as erroneous unaudited proposals divided by all proposals. These denominators answer different questions. Mark the first undefined when no proposal errors exist. If audited corrections can fail, also count erroneous final decisions, including audited ones.

   With selective labels, mark realized escaped-error counts and population error rates unknown. Optional inverse-probability estimates require valid logged inclusion probabilities and completed feedback; identify them as estimates and state their uncertainty method. Never present a ratio of estimated totals as an observed rate or an automatically certified bound.

6. **Validate and document evidence.** Replay a small fully labeled sequence while exposing labels to policy updates only on selected audits. Check that high-confidence cases retain spot checks, an unaudited round leaves the basic threshold unchanged, corrections preserve original labels, and rerunning logged draws reproduces routing. Compare workload and errors on the same sequence; separate replay results from live observations.

   The paper's experiments reuse fixed streams across randomized trials; qualitative annotation uses a model auditor proxy, and some results exceed target error levels. These limitations are reported in [Section 6](https://arxiv.org/html/2610.31662v1#S6).

   Produce a guarantee ledger listing the exact source version and theorem, controlled metric, finite-sample term, feedback assumptions, reviewer reliability, randomness, implementation deviations, and evidence status. Do not certify human correctness from reviewer agreement alone. State whether the output is a policy specification, a validated replay, or an operational monitoring report.

## Inputs / Outputs

**Inputs:** Ordered proposals and confidence scores; model and verifier versions; an error rubric; user-selected `alpha`; policy parameters; expert feedback or a labeled replay set; review timing and capacity constraints; and authorized output destination. Request missing essentials before claiming a working monitor, while preparing the schema and evidence ledger independently.

**Outputs:** A monitoring contract, versioned routing policy, decision/audit ledger, metric report with explicit denominators and unknowns, and guarantee ledger with source links. A valid 15-minute execution produces these artifacts for a bounded sample or reports concrete missing inputs; it does not establish long-run convergence within that time.

**Completion check:** Every routed sample has a score or explicit exception, a recorded policy and audit probability, a review status, and a traceable final disposition. Report metrics only where their labels and denominators support calculation. Keep provisional or pending decisions visible, and distinguish empirical observations from mathematical guarantees.

## References

- [Reliable AI-Assisted Decisions via Selective Auditing, arXiv:2610.31662v1](http://arxiv.org/abs/2610.31662v1). Isaac Gibbs, Elise L. Han, Ryan J. Tibshirani, and Tijana Zrnic; submitted October 22, 2026. Candidate source checked against the [versioned full text](https://arxiv.org/html/2610.31662v1), particularly Sections 2–3, 5–7. Operational logging and completion checks above are workflow recommendations, not additional theorems from the paper.
