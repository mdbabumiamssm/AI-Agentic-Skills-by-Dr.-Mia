<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Multi-Agent Systems (MAS)

**Category:** Agentic AI / Architectures
**Difficulty:** Intermediate/Advanced

## Overview
Single agents are powerful, but **Multi-Agent Systems** allow for specialized roles, better error checking, and "wisdom of the crowd." 

## Common Design Patterns

### 1. Supervisor / Coordinator
A central agent (The Supervisor) breaks down a complex task and delegates sub-tasks to specialized worker agents (e.g., "Researcher", "Coder", "Reviewer"). The Supervisor aggregates the results.

### 2. Debate
Two or more agents with conflicting personas (e.g., "Proponent" vs. "Opponent") argue a point. This forces the model to hallucinate less and uncover hidden nuances. A "Judge" agent (or the user) evaluates the final arguments.

### 3. Swarm
A large number of simple agents follow simple local rules to achieve complex emergent behavior, similar to ants or bees.

## Implementation (`debate_supervisor.py`)
This script implements a **Debate** pattern managed by a **Supervisor**.
-   **Optimist Agent:** Always looks for benefits.
-   **Critic Agent:** Always looks for flaws.
-   **Supervisor:** Manages the flow of conversation and summarizes the result.

## Libraries
In production, use:
-   **LangGraph:** For stateful, graph-based agent orchestration.
-   **AutoGen:** Microsoft's framework for conversable agents.
-   **CrewAI:** For role-playing agent teams.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->