# Agentic AI Resource Map (April 16, 2026)

This document is the high-signal source map for keeping the repository current. It is intentionally opinionated: official docs and official repositories come first, and benchmark papers are included only where they still help explain why a pattern matters.

## First-tier framework references

Use these as the default research starting points for greenfield agent builds:

- OpenAI building agents guide: https://platform.openai.com/docs/guides/agents
- OpenAI Agents SDK guide: https://platform.openai.com/docs/guides/agents-sdk
- OpenAI Agents SDK repo: https://github.com/openai/openai-agents-python
- LangGraph overview: https://docs.langchain.com/oss/python/langgraph/overview
- LangGraph repo: https://github.com/langchain-ai/langgraph
- PydanticAI docs: https://ai.pydantic.dev/
- PydanticAI repo: https://github.com/pydantic/pydantic-ai
- Google ADK docs: https://google.github.io/adk-docs/
- Google ADK repo: https://github.com/google/adk-python
- Microsoft Agent Framework docs: https://learn.microsoft.com/en-us/agent-framework/get-started/
- Microsoft Agent Framework repo: https://github.com/microsoft/agent-framework

## Strong secondary framework references

These are active and worth learning, but should be chosen for fit rather than habit:

- Agno docs: https://docs.agno.com/
- Agno repo: https://github.com/agno-agi/agno
- CrewAI docs: https://docs.crewai.com/en
- CrewAI repo: https://github.com/crewAIInc/crewAI
- Mastra docs: https://mastra.ai/docs
- Mastra repo: https://github.com/mastra-ai/mastra
- LlamaIndex Workflows docs: https://developers.llamaindex.ai/python/llamaagents/workflows/
- LlamaIndex repo: https://github.com/run-llama/llama_index
- smolagents docs: https://huggingface.co/docs/smolagents/en/index
- smolagents repo: https://github.com/huggingface/smolagents
- DSPy docs: https://dspy.ai/
- DSPy repo: https://github.com/stanfordnlp/dspy

## Evaluation and observability anchors

These should sit next to framework docs because they determine whether a workflow is deployable:

- OpenAI evals guide: https://platform.openai.com/docs/guides/evals
- OpenAI Agents SDK tracing docs: https://openai.github.io/openai-agents-python/tracing/
- LangSmith observability docs: https://docs.langchain.com/langsmith/observability
- LangSmith evaluation docs: https://docs.langchain.com/langsmith/evaluation
- Pydantic Logfire docs: https://logfire.pydantic.dev/
- CrewAI observability docs: https://docs.crewai.com/en/observability

## Specialized browser and local execution references

Use these when real browser or local machine execution is central to the workload rather than incidental:

- Browser Use docs: https://docs.browser-use.com/cloud/quickstart
- Browser Use workspaces docs: https://docs.browser-use.com/cloud/agent/workspaces
- Browser Use repo: https://github.com/browser-use/browser-use
- OpenClaw site: https://openclaw.ai
- OpenClaw org: https://github.com/openclaw
- OpenClaw repo: https://github.com/openclaw/openclaw
- OSWorld benchmark: https://os-world.github.io/

## Interoperability standards

These should anchor repository guidance whenever tools and remote agents must remain portable:

- Model Context Protocol introduction: https://modelcontextprotocol.io/introduction
- Model Context Protocol specification: https://modelcontextprotocol.io/specification/2025-11-25/basic/index
- MCP servers repo: https://github.com/modelcontextprotocol/servers
- MCP registry repo: https://github.com/modelcontextprotocol/registry
- Agent2Agent protocol docs: https://a2aproject.github.io/A2A/latest/
- A2A Python SDK repo: https://github.com/a2aproject/a2a-python
- A2A JavaScript SDK repo: https://github.com/a2aproject/a2a-js

## Provider operations anchors

Use official provider docs and SDK repos before blog summaries or community wrappers:

- OpenAI docs: https://platform.openai.com/docs/
- OpenAI deprecations: https://platform.openai.com/docs/deprecations
- OpenAI Python SDK: https://github.com/openai/openai-python
- OpenAI JavaScript SDK: https://github.com/openai/openai-node
- Anthropic docs: https://docs.anthropic.com/
- Anthropic model overview: https://docs.anthropic.com/en/docs/about-claude/models/overview
- Anthropic Python SDK: https://github.com/anthropics/anthropic-sdk-python
- Google Gemini docs: https://ai.google.dev/gemini-api/docs
- Google GenAI Python SDK: https://github.com/googleapis/python-genai
- AWS Bedrock docs: https://docs.aws.amazon.com/bedrock/
- AWS Bedrock samples: https://github.com/aws-samples/amazon-bedrock-samples
- Azure AI Foundry docs: https://learn.microsoft.com/en-us/azure/ai-foundry/
- Azure SDK for Python: https://github.com/Azure/azure-sdk-for-python
- Cohere docs: https://docs.cohere.com/
- Cohere Python SDK: https://github.com/cohere-ai/cohere-python
- Mistral docs: https://docs.mistral.ai/
- Mistral Python SDK: https://github.com/mistralai/client-python
- DeepSeek docs: https://api-docs.deepseek.com/
- DeepSeek org: https://github.com/deepseek-ai
- xAI docs: https://docs.x.ai/
- xAI Python SDK: https://github.com/xai-org/xai-sdk-python
- xAI cookbook: https://github.com/xai-org/xai-cookbook

## Transition material references

These should not be treated as first-choice stacks, but they still matter for migration and historical context:

- OpenAI Swarm repo: https://github.com/openai/swarm
- OpenAI Agents SDK repo: https://github.com/openai/openai-agents-python

## Foundational and benchmark papers

Keep these close to the framework docs so implementation choices stay grounded in actual task patterns:

- ReAct: https://arxiv.org/abs/2210.03629
- Toolformer: https://arxiv.org/abs/2302.04761
- Reflexion: https://arxiv.org/abs/2303.11366
- Tree of Thoughts: https://arxiv.org/abs/2305.10601
- Voyager: https://arxiv.org/abs/2305.16291
- WebArena: https://arxiv.org/abs/2307.13854
- Language Agent Tree Search: https://arxiv.org/abs/2310.04406
- GAIA benchmark: https://arxiv.org/abs/2311.12983
- SWE-agent: https://arxiv.org/abs/2405.15793
- DSPy paper: https://arxiv.org/abs/2310.03714
- OSWorld paper and benchmark: https://os-world.github.io/

## Repository curation rules

- Prefer one strong skill with a clear selection boundary over three overlapping skills.
- Prefer official docs and official repos for current operational facts.
- Add benchmark papers when they clarify why a framework pattern matters, not to make the repo look academic.
- Mark older but still-active frameworks as transition material when they are no longer the best default.
- Refresh index files whenever a new first-tier, observability, or execution-surface skill is introduced.
