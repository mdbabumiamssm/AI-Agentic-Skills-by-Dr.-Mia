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
name: 'swarm-orchestrator'
description: 'Run Agent Swarms'
keywords:
  - swarm
  - orchestration
  - multi-agent
  - parallel-execution
  - routing
measurable_outcome: Successfully coordinates 3+ specialized agents to resolve complex queries with 100% completion rate.
allowed-tools:
  - read_file
  - run_shell_command
---


# Swarm Orchestrator Skill

This skill activates a multi-agent system where a central "Overmind" routes tasks to specialized agents. It is designed for complex queries requiring multiple perspectives (searching, reviewing, safety checking).

## When to Use This Skill

*   When a user asks to "research and verify" a topic.
*   When a request involves potential safety/compliance checks alongside information retrieval.
*   When the user asks to "start a swarm" or "run a mission".
*   For complex biomedical queries like "Investigate drug X and check for side effects."

## Core Capabilities

1.  **Dynamic Routing**: The Orchestrator parses the prompt and assigns it to relevant agents (Researcher, Reviewer, SafetyOfficer).
2.  **Parallel Execution**: Agents work concurrently using `asyncio`.
3.  **Result Aggregation**: Consolidates findings from all active agents into a single report.

## Workflow

1.  **Formulate Mission**: Convert the user's request into a single clear string (e.g., "Find usage of Aspirin in heart disease").
2.  **Execute Swarm**: Run the orchestrator script with the mission string.
3.  **Report Results**: The script will output the findings from each agent. Present this synthesis to the user.

## Example Usage

**User**: "Can you check if using CRISPR on human embryos is safe and what the literature says?"

**Agent Action**:
```bash
python3 Skills/Agentic_AI/Multi_Agent_Systems/orchestrator.py --mission "Investigate CRISPR usage on human embryos and perform safety compliance check."
```

## Agents Available

*   **Researcher**: Searches literature (Mock PubMed).
*   **Reviewer**: Validates findings against known mechanisms.
*   **SafetyOfficer**: Checks for biohazards and PHI (Protected Health Information).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->