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
name: 'multi-turn-coding-agent-evaluation'
description: 'Evaluate coding agents with SWE-DIVE dynamic user simulation, conversation validation, retained transcripts, and final code tests under matched interaction and resource budgets.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Multi-Turn Coding Agent Evaluation

## Overview

Evaluate coding agents in sessions where users clarify requests and introduce requirements across multiple turns. SWE-DIVE combines a dynamic user simulator, a conversation validator, and Harbor orchestration with executable task tests; this skill turns that architecture into a reproducible evaluation workflow. Treat results as benchmark evidence, with explicit limits on claims about production performance.

## When to Use This Skill

- Reproduce a published SWE-DIVE configuration or run a clearly labeled adaptation.
- Compare coding agents on evolving requirements under matched budgets.
- Diagnose whether failure arose from code changes, user simulation, conversation validity, or infrastructure.
- Audit interaction transcripts alongside final code correctness.

## Core Capabilities

1. **Establish a reproducible configuration.** Read the versioned paper and repository documentation linked in References. Identify the published task selection, simulator behaviors, validator procedure, and Harbor integration before executing. Pin the repository commit, task revisions, environment image or dependency lock, model identifiers, prompts, and sampling settings. Use commands and configuration keys verified in that checkout; do not invent a launcher or claim an exact reproduction when required artifacts are unavailable.

2. **Define a bounded evaluation.** Select tasks and repetitions before observing outcomes. Record per-episode turn, token, tool-call, wall-clock, and cost limits where measurable, including separate simulator and validator consumption. For a 15-minute invocation, use an already prepared environment and a small smoke evaluation, reserving time for final tests and reporting. A full benchmark may require a separate longer run; record a blocked or incomplete outcome when setup or execution exceeds the window, without presenting it as completed evaluation.

3. **Execute dynamic interactions.** Start each episode from the same clean task revision and initial requirements. Let the simulator respond to the actual agent history using the selected behavior configuration, retaining each user message, agent response, tool event, and termination reason in order. Record requirement additions and clarifications as they occur. Keep evaluation-only tests and reference solutions outside agent-visible context according to the published protocol.

4. **Validate the simulated conversation.** Apply the published validator at its documented checkpoints and retain its decisions, rationale, and any rejected simulator outputs. Check whether simulator messages remain consistent with the task requirements and guide a solvable interaction. Follow the documented handling of invalid conversations; do not silently regenerate until success. If a custom validator or retry policy is necessary, label it as an adaptation and record its budget and effect on eligibility.

5. **Check final code independently.** Preserve the final patch, then run the benchmark's authoritative test procedure in the pinned environment. Capture commands, exit codes, test results, and execution errors. Keep conversation validity and test correctness as separate fields: validator approval does not establish code correctness, and passing available tests does not establish coverage of every conversational requirement. Map requirements to test evidence or mark them unverified.

6. **Compare agents fairly.** Match task IDs, initial revisions, simulator and validator versions, behavior profiles, seeds where supported, environment access, stopping rules, and resource ceilings. Adaptive conversations may diverge; hold the simulation policy fixed rather than claiming identical dialogue. Record actual usage as well as ceilings, and account for retries consistently. Report paired task outcomes and repetition counts; distinguish a one-shot baseline, if requested, as a different interaction condition.

7. **Produce an auditable result.** Report all attempted episodes, valid conversations, test successes and failures, invalid simulations, timeouts, and infrastructure failures. State the denominator for every rate and expose exclusions. Preserve the benchmark's official scoring when available; label additional summaries as local analysis. Explain simulator dependence, task coverage, stochastic variation, and any deviations from the publication without inferring production reliability.

These operating controls are recommended evaluation practices; only settings verified against the paper and pinned implementation should be described as the published SWE-DIVE protocol.

## Inputs / Outputs

**Inputs**

- Paper version and SWE-DIVE repository commit, with a runnable Harbor environment.
- Task IDs, initial repository revisions, initial requirements, and authoritative test specifications.
- Agent configurations and available credentials supplied through the environment.
- Simulator behavior profiles, validator settings, prompts, seeds, and repetition policy.
- Interaction and resource limits, requested comparison conditions, and an artifact destination.

**Outputs**

- A run manifest containing configuration versions, exact commands, budgets, deviations, and start/end times.
- Ordered per-episode transcripts with role, turn index, timestamp, message or tool event, and available usage counts; redact credentials while retaining evaluation-relevant content.
- Validator decisions, requirement-to-evidence mappings, final patches, and test logs linked by episode ID.
- A result table with at least `episode_id`, `task_id`, `agent_id`, `behavior_profile`, `seed`, `conversation_validity`, `test_status`, `termination_reason`, `turns`, `usage`, and artifact locations.
- A comparison report with explicit counts and denominators, observed resource use, unresolved evidence gaps, and completion status.

**Completion check:** A completed smoke evaluation must link its manifest, transcript, validator outcome, final patch, and executed test evidence. Failed agent tests are valid evaluation results. Missing prerequisites, unexecuted tests, or exhausted budgets require an explicit blocked or incomplete status; configuration preparation alone is not a successful reproduction.

## References

- [SWE-DIVE: Multi-Turn Coding Agent Evaluation with Dynamic User Simulation, arXiv:2610.30047v1](http://arxiv.org/abs/2610.30047v1) — Lorenz Brehme, Benjamin Steenhoek, and Ruth Breu; source finding, submitted October 21, 2026.
- [SWE-DIVE source repository](https://github.com/lorenzbrehme/swe-dive_usersim) — implementation and accompanying artifacts; inspect the pinned revision for executable setup and evaluation instructions.
