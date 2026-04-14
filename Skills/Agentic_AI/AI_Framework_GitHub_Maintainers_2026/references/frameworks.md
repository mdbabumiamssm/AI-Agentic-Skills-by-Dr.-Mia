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

# Official Framework and Protocol References (checked April 13, 2026)

## Recommended First-Tier Agent Frameworks

- OpenAI Agents SDK docs: https://platform.openai.com/docs/guides/agents-sdk
- OpenAI Agents SDK repo: https://github.com/openai/openai-agents-python
- LangGraph docs: https://docs.langchain.com/oss/python/langgraph/overview
- LangGraph repo: https://github.com/langchain-ai/langgraph
- LangSmith observability docs: https://docs.langchain.com/langsmith/observability
- PydanticAI docs: https://ai.pydantic.dev/
- PydanticAI repo: https://github.com/pydantic/pydantic-ai
- Google ADK docs: https://google.github.io/adk-docs/
- Google ADK Python repo: https://github.com/google/adk-python
- Microsoft Agent Framework docs: https://learn.microsoft.com/en-us/agent-framework/get-started/
- Microsoft Agent Framework repo: https://github.com/microsoft/agent-framework

## Strong Secondary Frameworks

- Agno docs: https://docs.agno.com/
- Agno repo: https://github.com/agno-agi/agno
- CrewAI docs: https://docs.crewai.com/
- CrewAI repo: https://github.com/crewAIInc/crewAI
- Mastra docs: https://mastra.ai/docs
- Mastra repo: https://github.com/mastra-ai/mastra
- LlamaIndex Workflows docs: https://developers.llamaindex.ai/python/llamaagents/workflows/
- LlamaIndex repo: https://github.com/run-llama/llama_index
- smolagents docs: https://huggingface.co/docs/smolagents
- smolagents repo: https://github.com/huggingface/smolagents
- DSPy docs: https://dspy.ai/
- DSPy repo: https://github.com/stanfordnlp/dspy

## Specialist Layers and Browser-Native Execution

- Browser Use quickstart: https://docs.browser-use.com/cloud/quickstart
- Browser Use workspaces docs: https://docs.browser-use.com/cloud/agent/workspaces
- Browser Use repo: https://github.com/browser-use/browser-use
- OSWorld benchmark: https://os-world.github.io/

## Interoperability Protocols

- Agent2Agent protocol docs: https://a2aproject.github.io/A2A/latest/
- Agent2Agent GitHub organization: https://github.com/a2aproject
- A2A Python SDK repo: https://github.com/a2aproject/a2a-python
- A2A JavaScript SDK repo: https://github.com/a2aproject/a2a-js
- Model Context Protocol overview: https://modelcontextprotocol.io/introduction
- Model Context Protocol specification: https://modelcontextprotocol.io/specification/2025-11-25/basic/index
- MCP servers repo: https://github.com/modelcontextprotocol/servers
- MCP registry repo: https://github.com/modelcontextprotocol/registry

## Transition and Legacy Notes

- AutoGen docs: https://microsoft.github.io/autogen/
- AutoGen repo: https://github.com/microsoft/autogen
- Semantic Kernel docs: https://learn.microsoft.com/en-us/semantic-kernel/
- Semantic Kernel repo: https://github.com/microsoft/semantic-kernel

### Current guidance

- For new Microsoft-centered builds, evaluate Microsoft Agent Framework before starting from AutoGen or Semantic Kernel.
- For durable long-running agents, evaluate LangGraph before assembling ad hoc orchestration around prompt loops.
- For strongly typed Python services and schema-heavy workflows, evaluate PydanticAI early.
- For Python teams that want one runtime spanning agents, teams, workflows, knowledge, and production serving, evaluate Agno when first-tier stacks do not fit.
- For role-based multi-agent automation with explicit crew and flow concepts, evaluate CrewAI as a strong secondary choice rather than a blind default.
- For TypeScript-first product teams that need agents, workflows, MCP, evals, and observability together, evaluate Mastra explicitly.
- For modular LM programming, prompt optimization, and evaluation-first improvement loops, evaluate DSPy as a specialist layer rather than a full runtime replacement.
- For real browser execution, keep Browser Use as the browser substrate and pair it with a larger orchestration framework when needed.
- For cross-agent interoperability between independently hosted agents, use A2A rather than inventing a custom RPC layer.
- For tool and resource interoperability, use MCP rather than framework-specific tool adapters when cross-client reuse matters.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
