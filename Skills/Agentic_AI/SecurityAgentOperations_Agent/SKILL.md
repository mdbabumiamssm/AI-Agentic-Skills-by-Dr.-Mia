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
name: 'security-agent-operations'
description: 'Operate security-focused coding agents with scoped authorization, lab isolation, adversarial workflows, tool boundaries, and defensive validation inspired by Raptor.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Security Agent Operations

## Overview

This skill guides safe operation of AI agents configured for security research, red-team analysis, blue-team validation, and controlled adversarial reasoning. It emphasizes explicit authorization, bounded scope, lab isolation, tool governance, and defensive evidence so agentic workflows remain useful without drifting into unauthorized or unsafe activity.

## When to Use This Skill

- Building or operating an AI security agent with rules, sub-agents, tools, skills, or task-specific playbooks.
- Planning vulnerability research, threat modeling, exploitability analysis, detection engineering, or incident-response support.
- Running offensive-style tests in an owned, contracted, CTF, sandbox, or lab environment.
- Converting adversarial findings into defensive validation, hardening steps, detections, or regression tests.
- Reviewing whether a proposed agent workflow has adequate legal scope, isolation, logging, and tool boundaries.

## Core Capabilities

1. Authorization and scope control: Confirm the asset owner, written permission, test window, target inventory, prohibited actions, data handling rules, and stop conditions before any active security work.

2. Lab isolation: Prefer disposable infrastructure, local fixtures, intentionally vulnerable targets, CTF ranges, or cloned test systems; avoid touching third-party production assets unless the authorization explicitly covers them.

3. Agent rule design: Define system rules that constrain objectives, forbid unauthorized access, require evidence capture, and force escalation when scope or legality is uncertain.

4. Tool boundary management: Classify tools by risk, require dry-run or read-only modes where possible, and gate intrusive scanners, exploit frameworks, credential tooling, traffic generators, and destructive commands behind explicit scope checks.

5. Adversarial workflow planning: Decompose work into reconnaissance, hypothesis generation, validation, impact assessment, remediation, and retest phases while keeping each phase tied to permitted assets.

6. Defensive validation: Translate findings into patch guidance, configuration changes, detection logic, log queries, test cases, and reproducible verification steps.

7. Evidence and auditability: Keep concise records of commands, timestamps, target identifiers, tool versions, observed outputs, screenshots when relevant, and the reasoning that links evidence to conclusions.

8. Safety stop conditions: Pause when encountering out-of-scope systems, live secrets, personal data, destructive side effects, unclear authorization, persistence mechanisms, or requests to evade monitoring.

9. Raptor-style security-agent operation: Frame adversarial tasks only inside scoped engagement files and isolated labs, separate red-team and blue-team task modes, bind sub-agents and tools to explicit allowlists, enforce adversarial reasoning rules, maintain evidence logs of actions and reasoning, and require defensive validation before reporting results or taking real-world security action.

## Inputs / Outputs

Inputs:

- Authorization statement, rules of engagement, target list, and testing window.
- Agent configuration files, prompt rules, tool manifests, sub-agent definitions, or security playbooks.
- Repository, system, container, VM, log, alert, or network-lab context relevant to the assessment.
- User objective such as threat model, vulnerability triage, defensive test plan, remediation validation, or incident-support workflow.

Outputs:

- Scope-checked operational plan with permitted assets, excluded actions, assumptions, and stop conditions.
- Agent rule recommendations, tool restrictions, and workflow checkpoints.
- Evidence-backed findings with severity rationale, reproduction steps limited to authorized environments, and defensive impact.
- Remediation, detection, regression-test, and retest guidance.
- Audit notes sufficient for handoff to engineering, security operations, or governance stakeholders.

## References

- Raptor: https://github.com/gadievron/raptor
- MITRE Caldera: https://github.com/mitre/caldera
- Atomic Red Team: https://github.com/redcanaryco/atomic-red-team
- OWASP Application Security Verification Standard: https://github.com/OWASP/ASVS
