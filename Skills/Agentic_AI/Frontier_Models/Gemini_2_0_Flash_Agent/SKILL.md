---
name: Gemini_2_0_Multimodal_Agent
description: High-speed multimodal agent capable of processing video, audio, and large context windows instantly.
version: 2.0.0
category: Agentic_AI
author: AI Agentic Skills Team
input:
  query: Prompt or question.
  media_path: Path to video/audio file (optional).
output:
  response: Multimodal analysis result.
---

# Gemini 2.0 Multimodal Agent

Built on **Gemini 2.0 Flash**, this agent is designed for speed and "native multimodal" understanding. It doesn't just see frames; it understands temporal dynamics in video and nuance in audio natively.

## Capabilities
1.  **1M+ Context Window:** Can ingest entire codebases or long videos.
2.  **Native Audio/Video:** Process meeting recordings or lab videos in real-time.
3.  **Sub-second Latency:** Optimized for interactive agentic loops.

## Usage
```python
from Skills.Agentic_AI.Frontier_Models.Gemini_2_0_Flash_Agent.gemini_agent import GeminiAgent

agent = GeminiAgent()
result = agent.analyze_video(
    video_path="./lab_experiment.mp4",
    prompt="Identify the exact timestamp where the chemical reaction changes color."
)
```
