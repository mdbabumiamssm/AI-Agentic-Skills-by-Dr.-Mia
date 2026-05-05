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
name: 'microsoft-prior-auth-multi-agent-accelerator'
description: 'Guide payer-side prior authorization review using Microsoft Azure multi-agent patterns, Foundry-hosted agents, MCP healthcare data access, audit trails, and human oversight.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Microsoft Prior Auth Multi-Agent Accelerator

## Overview

Use this skill to evaluate, adapt, or explain the Microsoft Prior Authorization Multi-Agent Solution Accelerator for payer-side AI-assisted prior authorization review. It focuses on the accelerator's multi-agent architecture, gate-based decision rubric, healthcare MCP data access, confidence scoring, auditability, and human-in-the-loop review patterns for enterprise clinical operations.

## When to Use This Skill

- Designing a payer-side prior authorization workflow using Microsoft Agent Framework or Azure AI Foundry-hosted agents.
- Comparing multi-agent prior authorization patterns against existing utilization management or clinical review operations.
- Mapping compliance, clinical, coverage, and synthesis responsibilities into separate review agents.
- Planning an Azure Container Apps deployment using `azd` for a prior authorization accelerator.
- Creating governance, audit trail, confidence scoring, or human oversight requirements for AI-assisted authorization review.
- Explaining how healthcare MCP access can support evidence retrieval during authorization review.

## Core Capabilities

1. **Architecture orientation** - Identify the accelerator's four-agent pattern: Compliance, Clinical, Coverage, and Synthesis agents coordinated for payer-side review.
2. **Review rubric mapping** - Translate gate-based decision logic into explicit review checkpoints for compliance, clinical necessity, benefit coverage, and final synthesis.
3. **Healthcare data access planning** - Describe how MCP-based healthcare data access can be used to retrieve records, coverage details, and policy context needed by review agents.
4. **Confidence and escalation design** - Use confidence scoring and uncertainty thresholds to decide when a case should move to human review.
5. **Audit trail requirements** - Capture inputs, retrieved evidence, agent outputs, rubric outcomes, confidence scores, and final reviewer actions for traceability.
6. **Azure deployment fit check** - Assess whether Azure Container Apps, Azure Developer CLI deployment, and Foundry-hosted agents fit the target enterprise environment.
7. **Clinical operations adaptation** - Convert the reference accelerator into an implementation plan that respects payer policies, clinical review workflows, and oversight responsibilities.

## Inputs / Outputs

**Inputs**

- Prior authorization case context, including requested service, diagnosis, patient history, payer policy, and coverage information.
- Existing clinical operations workflow, escalation rules, and human reviewer responsibilities.
- Target Azure environment assumptions, including Azure AI Foundry, Azure Container Apps, and Azure Developer CLI readiness.
- Compliance, audit, privacy, and governance requirements for AI-assisted healthcare review.

**Outputs**

- A concise architecture summary of the Microsoft multi-agent prior authorization pattern.
- A case-review workflow assigning responsibilities to Compliance, Clinical, Coverage, and Synthesis agents.
- A gate-based decision rubric with required evidence, confidence handling, and human escalation criteria.
- An implementation or evaluation checklist for Azure deployment and clinical operations adoption.
- An audit trail schema or checklist covering evidence, agent reasoning summaries, decisions, confidence scores, and reviewer actions.

## References

- Microsoft Prior Authorization Multi-Agent Solution Accelerator: https://github.com/microsoft/Prior-Authorization-Multi-Agent-Solution-Accelerator
