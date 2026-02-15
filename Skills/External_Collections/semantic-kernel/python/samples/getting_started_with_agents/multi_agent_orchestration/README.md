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

# Multi-agent orchestration

The Semantic Kernel Agent Framework now supports orchestrating multiple agents to work together to complete a task.

## Background

The following samples are beneficial if you are just getting started with Semantic Kernel.

- [Chat Completion](../../concepts/chat_completion/)
- [Auto Function Calling](../../concepts/auto_function_calling/)
- [Structured Output](../../concepts/structured_output/)
- [Getting Started with Agents](../../getting_started_with_agents/)
- [More advanced agent samples](../../concepts/agents/)

## Prerequisites

The following environment variables are required to run the samples:

- OPENAI_API_KEY
- OPENAI_CHAT_MODEL_ID

However, if you are using other model services, feel free to switch to those in the samples.
Refer to [here](../../concepts/setup/README.md) on how to set up the environment variables for your model service.

## Orchestrations

| **Orchestrations** | **Description**                                                                                                                                                                                     |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Concurrent**     | Useful for tasks that will benefit from independent analysis from multiple agents.                                                                                                                  |
| **Sequential**     | Useful for tasks that require a well-defined step-by-step approach.                                                                                                                                 |
| **Handoff**        | Useful for tasks that are dynamic in nature and don't have a well-defined step-by-step approach.                                                                                                    |
| **GroupChat**      | Useful for tasks that will benefit from inputs from multiple agents and a highly configurable conversation flow.                                                                                    |
| **Magentic**   | GroupChat like with a planner based manager. Inspired by [Magentic One](https://www.microsoft.com/en-us/research/articles/magentic-one-a-generalist-multi-agent-system-for-solving-complex-tasks/). |


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->