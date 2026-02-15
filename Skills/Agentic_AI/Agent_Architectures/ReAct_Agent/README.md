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

# ReAct Agent (Reasoning + Acting)

**Category:** Agentic AI / Architectures
**Difficulty:** Intermediate

## Overview
The **ReAct (Reasoning + Acting)** pattern is a foundational paradigm in Agentic AI. Instead of just answering a question, a ReAct agent explicitly generates a "Thought" (reasoning trace), decides on an "Action" (tool use), receives an "Observation" (tool output), and repeats this cycle until it arrives at a "Final Answer".

Paper: [ReAct: Synergizing Reasoning and Acting in Language Models](https://arxiv.org/abs/2210.03629)

## Components
1.  **Agent Loop:** The `think()` method in `react_core.py` orchestrates the Thought-Action-Observation loop.
2.  **Tools:** A dictionary of callable functions (defined in `tools.py`).
3.  **Prompting:** The agent requires a specific system prompt to enforce the `Thought -> Action -> Observation` format.

## Implementation Details

The core loop works as follows:
1.  Append user input to history.
2.  Send history to LLM.
3.  LLM generates: `Thought: ... Action: ToolName(Args)`
4.  System parses `Action`.
5.  System executes tool and gets `Observation`.
6.  `Observation` is appended to history.
7.  Repeat until `Final Answer` is generated.

## Usage

```python
from react_core import Agent
from tools import registry

# Mock LLM Client for demonstration
class MockLLM:
    def generate(self, prompt):
        # Simple hardcoded behavior for testing
        if "London" in prompt and "Observation" not in prompt:
            return "Thought: I need to check the weather in London. Action: Weather(London)"
        if "Observation: Rainy" in prompt:
            return "Thought: The weather is rainy. Final Answer: It is rainy in London."
        return "Final Answer: I don't know."

agent = Agent(MockLLM(), registry)
response = agent.think("What is the weather in London?")
print(response)
```

## Future Improvements
- Add support for real LLM APIs (OpenAI, Anthropic).
- Implement a more robust parser (handle parsing errors gracefully).
- Add "Memory" to persist context across sessions.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->