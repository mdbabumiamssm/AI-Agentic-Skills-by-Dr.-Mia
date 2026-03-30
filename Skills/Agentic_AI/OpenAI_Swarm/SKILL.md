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

---
name: 'openai-swarm-orchestrator'
description: 'Lightweight agent swarms built on OpenAI\'s Swarm framework (agents, handoffs, routines).'
keywords:
  - openai
  - swarm
  - handoffs
  - routines
  - lightweight-multi-agent
measurable_outcome: Enables agent teams with orchestrated handoffs that solve end-to-end biomedical investigations with <60 s orchestration overhead.
allowed-tools:
  - read_file
  - run_shell_command
  - python
---

# OpenAI Swarm Skill

OpenAI Swarm is a minimalist multi-agent framework centered on **explicit handoffs** between constrained, tool-calling agents. Agents encapsulate instructions + Python tool lists, and any tool can hand off execution by returning another `Agent` definition. This keeps multi-agent reasoning transparent and debuggable versus monolithic prompts.

* Core primitives ([openai/swarm](https://github.com/openai/swarm) overview, SourcePulse analysis). [^sourcepulse]
* Production guardrails & deploy patterns (Lexogrine 2026 Swarm field guide). [^lexogrine]

## When to use this skill

* You need rapid experiments with handoffs before migrating to the OpenAI Agents SDK.
* Regulatory workflows that demand deterministic routing + audit logs.
* Tool-heavy biomedical assistants where each expert agent maps 1:1 to a validated SOP.

## Quickstart

```bash
pip install --upgrade openai swarm
export OPENAI_API_KEY=sk-...
```

```python
from swarm import Agent, Swarm

triage = Agent(
    name="ClinicalTriage",
    instructions="Collect vitals + labs, produce triage_summary JSON.",
)

safety = Agent(
    name="SafetyOfficer",
    instructions="Check outputs against boxed warnings, escalate if needed.",
)

def handoff_to_safety(state):
    return safety  # explicit handoff once triage finds red flags

triage.functions = [handoff_to_safety]
result = Swarm().run(agent=triage, messages=[{"role": "user", "content": "Patient chest pain"}])
```

## Operational patterns

1. **Routines for expensive loops** – reuse `swarm.routines` to schedule repeated lab lookups.
2. **Context injection** – pass `context_variables` so downstream agents know patient_id/mission_id.
3. **Eval bridge** – wrap Swarm runs with `platform/evaluators/swarm_eval.py` to verify SOP compliance.
4. **Upgrade path** – once workflows stabilize, migrate to Agents SDK (see `AI_Providers/OpenAI_Agents_SDK_2026`).

## Deliverables

* `playbooks/swarm_blueprints.md` – canonical agent templates (triage, safety, automation) for internal review.
* `scripts/run_swarm_demo.py` – runnable example connecting Swarm to in-repo bioscience tools.

## References

[^sourcepulse]: SourcePulse project summary for [openai/swarm](https://www.sourcepulse.org/projects/1832643) highlighting agent primitives and handoffs.
[^lexogrine]: Lexogrine (Feb 16, 2026) "OpenAI Swarm Multi-Agent Framework in 2026: What It Is, How It Works, and How to Use It" explaining handoffs, approval gates, and upgrade guidance.

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
