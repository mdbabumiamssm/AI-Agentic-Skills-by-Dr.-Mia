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

# Tutorial: Building a Bio-Medical Agent Swarm

**Level:** Intermediate
**Prerequisites:** Python 3.9+, AsyncIO basics
**Estimated Time:** 20 Minutes

## Introduction

Traditional automation relies on linear scripts. **Agentic Swarms** rely on specialized "workers" that can run in parallel, make decisions, and communicate. 

In this tutorial, we will use the `SwarmOrchestrator` to build a team of agents that can:
1.  **Search** medical literature.
2.  **Review** findings for accuracy.
3.  **Check** for safety/compliance.

## Core Concepts

### 1. The Orchestrator (The Brain)
The Orchestrator doesn't do the work. It **routes** the work. It takes a high-level objective (e.g., "Analyze this drug") and decides *who* needs to be involved.

### 2. The Agents (The Workers)
Each agent has a `role` and `capabilities`.
*   **ResearchAgent:** specialized in fetching data.
*   **ReviewAgent:** specialized in verifying data.
*   **SafetyAgent:** specialized in compliance.

### 3. Async Parallelism
Unlike a waterfall process (Step 1 -> Step 2), a Swarm can fire multiple agents at once. The `asyncio.gather` pattern allows the Researcher and Safety Officer to work simultaneously.

## Step-by-Step Implementation

### Step 1: Define your Agents
Inherit from `BaseAgent` and implement the `process` method.

```python
class MySpecialAgent(BaseAgent):
    async def process(self, message: AgentMessage) -> AgentMessage:
        # 1. Think (Call LLM)
        # 2. Act (Call Tool)
        return AgentMessage(...)
```

### Step 2: Register with the Swarm
```python
swarm = SwarmOrchestrator()
swarm.register_agent(MySpecialAgent("Agent007", "Spy", ["sneak"]))
```

### Step 3: Run a Mission
```python
await swarm.run_mission("Find the secret formula")
```

## Advanced: adding a Real LLM Router

The current `orchestrator.py` uses keyword matching (`if "search" in task`). To make this "Cognitive", replace `_route_task` with an LLM call:

```python
async def _route_task(self, task: Task) -> List[str]:
    system_prompt = f"You are a project manager. Available agents: {list(self.agents.keys())}."
    user_prompt = f"Which agents are needed for: {task.description}? Return JSON list."
    
    # Call OpenAI/Anthropic
    response = await call_llm(system_prompt, user_prompt)
    return json.loads(response) # e.g. ["Researcher", "SafetyOfficer"]
```

## Next Steps

1.  Open `orchestrator.py` and run it to see the output.
2.  Try adding a `ClinicalTrialAgent` that mocks checking ClinicalTrials.gov.
3.  Modify the "Mission" string to trigger your new agent.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->