---
name: DeepSeek_R1_Open_Agent
description: Open-weights reasoning model specialized in math, code, and logic with transparent Chain-of-Thought (CoT).
version: 1.0.0
category: Agentic_AI
author: AI Agentic Skills Team
input:
  query: Mathematical proof, complex algorithm, or logic puzzle.
output:
  reasoning: The raw Chain-of-Thought output.
  final_answer: The concise conclusion.
---

# DeepSeek R1 Open Agent

This agent wraps the **DeepSeek R1** model, a frontier-class open model that rivals proprietary models in reasoning tasks. It is unique for providing its raw "internal monologue" (Chain-of-Thought) to the user.

## Capabilities
1.  **Transparent Reasoning:** Returns the full `<think>` tags showing how it arrived at an answer.
2.  **Math & Logic:** Exceptional performance on AIME and Codeforces benchmarks.
3.  **Cost Efficiency:** Designed to be distilled into smaller models or run efficiently.

## Usage
```python
from Skills.Agentic_AI.Frontier_Models.DeepSeek_R1_Agent.deepseek_agent import DeepSeekAgent

agent = DeepSeekAgent()
result = agent.solve(
    problem="Prove that there are infinitely many prime numbers of the form 4n+3."
)
print(result["reasoning"])
print(result["answer"])
```
