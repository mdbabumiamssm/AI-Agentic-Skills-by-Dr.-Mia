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
name: 'collaborative-agent-reasoning-engineering'
description: 'Design domain AI agents through CARE stage gates with SMEs, developers, helper agents, reviewable specs, reasoning policies, tools, and evaluations.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Collaborative Agent Reasoning Engineering

## Overview

Use this skill to engineer domain AI agents with the Collaborative Agent Reasoning Engineering (CARE) methodology. CARE turns informal subject-matter intent into reviewable artifacts, stage-gated decisions, tool orchestration rules, reasoning policies, and evaluation criteria so agent behavior can be specified, tested, and maintained.

This skill is most useful when building agents for scientific or other constrained domains where expert approval, grounded behavior, reproducible workflows, and verification practices matter more than quick prompt iteration.

## When to Use This Skill

- A user wants to design, specify, or improve a domain AI agent with SME review.
- Agent behavior needs explicit interaction requirements, reasoning policies, grounding rules, or tool-use constraints.
- The work involves translating informal domain intent into structured, reviewable implementation artifacts.
- A team needs stage gates before implementation, deployment, or evaluation.
- The agent must operate in a scientific, technical, biomedical, compliance-sensitive, or high-verification workflow.
- Existing runtime framework guidance is not enough because the problem is collaborative specification and governance.

## Core Capabilities

1. **Stakeholder role mapping**: Identify the SME, developer, and helper-agent roles, then define who provides domain constraints, who implements the system, and who facilitates artifact generation.
2. **Intent capture**: Convert informal goals, workflows, data assumptions, and failure concerns into structured agent requirements for human review.
3. **Interaction requirements**: Specify user types, supported tasks, expected inputs, response boundaries, escalation paths, and unacceptable behaviors.
4. **Reasoning policy design**: Define how the agent should reason, cite evidence, handle uncertainty, ask clarifying questions, and avoid unsupported conclusions.
5. **Grounding and tool orchestration**: Document required sources, retrieval behavior, tool preconditions, tool sequencing, permission boundaries, and fallback behavior.
6. **Stage-gated artifact review**: Produce artifacts that SMEs and developers can approve before moving to implementation, testing, or deployment.
7. **Evaluation criteria**: Define test cases, acceptance checks, failure modes, review rubrics, and regression criteria tied to the agreed requirements.
8. **Maintainability planning**: Track assumptions, open questions, versioned artifacts, and change triggers so the agent can evolve without losing domain alignment.

## Inputs / Outputs

**Inputs**

- Domain goal, intended users, and target workflows.
- SME-provided constraints, terminology, examples, edge cases, and unacceptable outputs.
- Available data sources, APIs, tools, retrieval systems, or execution environments.
- Developer constraints such as runtime framework, deployment target, logging, privacy, or review process.
- Existing prompts, agent traces, evaluation results, or failure reports when improving an existing agent.

**Outputs**

- CARE role map describing SME, developer, and helper-agent responsibilities.
- Interaction requirements covering user intents, input assumptions, outputs, limits, and escalation behavior.
- Reasoning policy for evidence use, uncertainty handling, verification, and refusal or deferral behavior.
- Grounding and tool orchestration specification with source hierarchy, tool contracts, and fallback rules.
- Stage-gate checklist for SME/developer review before implementation, evaluation, and deployment.
- Evaluation plan with representative cases, expected behavior, failure checks, and maintenance triggers.

## References

- CARE source finding: http://arxiv.org/abs/2604.28043v1
