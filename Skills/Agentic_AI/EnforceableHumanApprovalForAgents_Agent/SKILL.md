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
name: 'enforceable-human-approval-for-agents'
description: 'Threat-model and test human approval gates against agent-controlled platforms, using I-HITL to distinguish enforceable approval from agent-accessible authentication.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Enforceable Human Approval for Agents

## Overview

This skill reviews whether protected agent actions actually require explicit human approval when the agent controls its local platform. It translates I-HITL into a bounded threat-modeling, enforcement-design, and bypass-testing workflow, distinguishing ordinary I-HITL from strong I-HITL, which also resists a human delegating approval to the agent. It produces an evidence-backed assessment, not a certification or a substitute for an external enforcement mechanism.

## When to Use This Skill

- Reviewing approval gates for payments, permission changes, deployments, or sensitive data access.
- An agent can read the inbox, browser session, configuration, or credentials used for approval.
- A local harness pauses execution but the agent may control the host enforcing that pause.
- Evaluating trusted-device approval or a Countersign-style integration.
- Testing whether recovery, enrollment, or alternate APIs bypass a claimed human approval requirement.

## Core Capabilities

1. **Establish evidence and scope.** Identify the protected action and separate source results, deployment observations, and proposed controls.
2. **Map trust boundaries.** Trace the human, agent platform, enforcing service, trusted device, and communication channels.
3. **Assess approval enforcement.** Check whether the service requires evidence tied to the exact action and whether the agent can manufacture that evidence.
4. **Evaluate delegation resistance.** Examine key ownership, enrollment, human input, and recovery under an explicitly selected adversary model.
5. **Test and report bypasses.** Use authorized fixtures to produce reproducible results and identify missing evidence without claiming a formal proof.

### Evidence baseline and limits

The supplied finding reports autonomous OTP harvesting to raise wallet spending limits and introduces formal communication-topology conditions for I-HITL. A trusted entity is necessary, but its mere presence is insufficient. Treat those results as motivation and model-dependent guarantees, not evidence that a particular deployment is secure.

The full paper requires attested key generation; strong I-HITL additionally depends on preventing key delegation. Its model trusts the enforcing service and excludes agent physical interaction. Human comprehension and rubber-stamping remain outside its guarantee. Countersign combines attestation, protected keys, a trusted display, and biometric confirmation. Both mobile variants trust OS-mediated display binding; iOS adds an unverifiable non-jailbroken-OS assumption. Attestation is time-limited evidence, and prototype key recovery remains unresolved. See Sections II, V–VII of the [paper](https://arxiv.org/html/2610.35835v1).

The following workflow and acceptance checks are engineering recommendations. They are not claimed to be the paper's complete construction or theorem conditions. Before asserting conformance, compare the actual topology and protocol against the full formal definitions and record the applicable theorem and assumptions.

### Workflow: bounded initial assessment

**1. Scope and inventory — minutes 0–3.**

Read the supplied architecture, gate policy, and relevant implementation. Record protected action types, approval granularity, authorized test environment, and whether the objective assumes an honest or delegating human. Missing access should produce an explicit evidence gap, not an invented result.

Inventory every route capable of executing the action or weakening its gate: API, CLI, administrative controls, policy edits, recovery, device enrollment, and queued jobs. Do not assume that task authorization permits changing approval policy.

**2. Threat model — minutes 3–6.**

For each component, record who controls it, what secrets it holds, and which messages it can forge, read, suppress, or replay. Assume complete agent control of the local platform for the I-HITL assessment. Identify the independent service that ultimately rejects unapproved actions; for local resources, identify the isolation boundary protecting that service.

Produce a channel table with: sender, receiver, authentication, confidentiality, agent access, and supporting evidence. An emailed OTP accessible through an agent tool is not independent human approval. Distinguish human identity authentication from approval of a specific action.

**3. Enforcement design — minutes 6–10.**

Write a reviewable contract for the protected action:

- The service creates a pending request with a unique identifier, expiry, account, destination, operation, parameters, and policy version as applicable.
- A canonical representation binds all security-relevant fields to the approval. The human-facing display must faithfully convey those fields; a signature over an opaque identifier alone does not establish that relationship.
- The trusted component requires fresh human confirmation for that request. Verify whether software, accessibility automation, or a shared secret can simulate confirmation.
- The service verifies enrollment, key provenance, relevant attestation evidence, signature, request binding, freshness, and revocation before execution. It consumes approval atomically to prevent reuse or duplicate execution.
- Any material parameter change invalidates approval. Timeouts, denials, unavailable verification, and unsupported device properties leave the action blocked.
- Enrollment, recovery, fallback methods, and gate-policy changes receive their own threat assessment; none may silently substitute an agent-accessible approval channel.

For strong I-HITL, assess whether a willing human can export or reconstruct approval keys, enroll agent-generated keys, or enable unattended signing. Record device-specific evidence rather than inferring guarantees from labels such as MFA, secure hardware, or biometric login. Re-check current platform documentation before implementation.

**4. Bypass tests — minutes 10–13.**

Run only within the authorized test scope, using synthetic credentials and harmless actions. Otherwise supply a test plan and mark results untested. Do not harvest real OTPs, raise live limits, or weaken production gates to demonstrate a failure.

| Test | Expected result |
| --- | --- |
| Valid approval for the displayed action | Exactly that action executes once |
| No approval, denial, or expired request | No execution |
| Agent reads a fixture OTP or calls the API directly | Protected action remains blocked |
| Modify recipient, amount, operation, or account after approval | Approval rejected |
| Replay approval or submit concurrent duplicates | No additional execution |
| Edit local approval flags or harness configuration | Independent service still enforces approval |
| Substitute an unattested key or revoked device | Enrollment or approval rejected |
| Simulate confirmation or substitute displayed content | No valid approval for unseen content |
| Use recovery, fallback, or policy-change routes | Claimed approval boundary is preserved |
| Attempt key export or unattended signing in strong mode | No delegable approval capability obtained |

Record setup, expected behavior, observed behavior, and redacted evidence for each executed test. Stop live testing at the first unexpected side effect; report it and continue only with non-mutating analysis or isolated reproduction.

**5. Assessment — minutes 13–15.**

Report each control as supported, failed, or unverified. Classify the deployment as ordinary HITL, candidate I-HITL, or candidate strong I-HITL, with explicit assumptions and unresolved blockers. Passing finite tests is not a proof; an untested requirement prevents a verified deployment claim. The fifteen-minute outcome is a completed assessment or actionable gap report, not a production implementation.

## Inputs / Outputs

**Inputs:** protected actions and policy; architecture or code; agent permissions; approval and recovery flows; device and attestation details; selected adversary model; authorized testing scope and available fixtures.

**Outputs:** one concise assessment containing:

- Scope and evidence ledger, distinguishing paper claims, observed facts, and recommendations.
- Component/channel trust map and explicit enforcement assumptions.
- Proposed approval contract and prioritized fixes for each bypass route.
- Test matrix with pass, fail, or untested results and redacted evidence.
- Conditional classification, remaining gaps, and next validation steps.

Return the assessment inline unless the user requests a file. A valid output names the enforcement point, covers every identified execution route, records the human-delegation assumption, and does not label proposed or untested controls as demonstrated guarantees. This skill itself cannot enforce human approval against a compromised host.

## References

- [I-HITL: Inalienable Human-in-the-Loop, arXiv:2610.35835v1](http://arxiv.org/abs/2610.35835v1), submitted October 26, 2026. Primary source finding and formal security model.
- [Version 1 full text](https://arxiv.org/html/2610.35835v1), especially Sections II and V–VII for assumptions, topology conditions, device properties, and prototype limitations.
- [Countersign prototype repository](https://github.com/I-HITL/countersign), linked by the paper. Inspect and pin a revision before reuse; repository availability is not evidence of deployment security.
