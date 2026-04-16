# Research Snapshot: Agentic AI Official Stack Health

Date: April 16, 2026

This file is a repository-local research artifact summarizing the upstream references used for the April 16, 2026 refresh.

## Official sources checked

### OpenAI

- Agents guide: https://platform.openai.com/docs/guides/agents
- Agents SDK guide: https://platform.openai.com/docs/guides/agents-sdk
- Evals guide: https://platform.openai.com/docs/guides/evals
- Agents SDK repo: https://github.com/openai/openai-agents-python
- Agents SDK releases: https://github.com/openai/openai-agents-python/releases
- Swarm repo: https://github.com/openai/swarm

Observed note:
- The official Agents SDK releases page shows `v0.14.1` on April 15, 2026.
- The official Swarm repo states that Swarm has been replaced by the OpenAI Agents SDK and recommends migrating production use cases.

### LangGraph and LangSmith

- LangGraph overview: https://docs.langchain.com/oss/python/langgraph/overview
- LangGraph repo: https://github.com/langchain-ai/langgraph
- LangGraph releases: https://github.com/langchain-ai/langgraph/releases
- LangSmith observability: https://docs.langchain.com/langsmith/observability
- LangSmith evaluation: https://docs.langchain.com/langsmith/evaluation

Observed note:
- The release stream is split across packages. `langgraph-checkpoint==4.0.2` appears on April 15, 2026, while alpha releases for the main package also appear in mid-April.

### PydanticAI and Logfire

- Docs: https://ai.pydantic.dev/
- Repo: https://github.com/pydantic/pydantic-ai
- Releases: https://github.com/pydantic/pydantic-ai/releases
- Logfire: https://logfire.pydantic.dev/

Observed note:
- The official releases page shows `v1.83.0` on April 16, 2026.

### Google ADK and Gemini

- ADK docs: https://google.github.io/adk-docs/
- ADK repo: https://github.com/google/adk-python
- ADK releases: https://github.com/google/adk-python/releases
- Gemini docs: https://ai.google.dev/gemini-api/docs
- Gemini models: https://ai.google.dev/gemini-api/docs/models
- Gemini release notes: https://ai.google.dev/gemini-api/docs/release-notes

Observed note:
- The official ADK releases page shows stable `v1.30.0` on April 13, 2026 and active `v2.0.0a3` pre-release work on April 9, 2026.

### Microsoft Agent Framework

- Get started: https://learn.microsoft.com/en-us/agent-framework/get-started/
- Repo: https://github.com/microsoft/agent-framework
- Releases: https://github.com/microsoft/agent-framework/releases

Observed note:
- The official releases page shows Python beta `1.0.0b260107` dated January 7, 2026.

### CrewAI and Agno

- CrewAI docs: https://docs.crewai.com/en
- CrewAI repo: https://github.com/crewAIInc/crewAI
- CrewAI releases: https://github.com/crewAIInc/crewAI/releases
- Agno docs: https://docs.agno.com/
- Agno repo: https://github.com/agno-agi/agno
- Agno releases: https://github.com/agno-agi/agno/releases

Observed note:
- CrewAI `1.8.0` and Agno `v2.3.24` appear on January 8, 2026 on their official releases pages.

### Browser Use and OpenClaw

- Browser Use docs: https://docs.browser-use.com/cloud/quickstart
- Browser Use repo: https://github.com/browser-use/browser-use
- Browser Use releases: https://github.com/browser-use/browser-use/releases
- OpenClaw org: https://github.com/openclaw
- OpenClaw repo: https://github.com/openclaw/openclaw
- OpenClaw releases: https://github.com/openclaw/openclaw/releases

Observed note:
- Browser Use `0.12.6` appears on April 2, 2026.
- Browser Use `0.12.5` documents a security-motivated dependency change after a March 2026 supply-chain incident affecting `litellm` versions.
- OpenClaw `2026.4.14` appears on April 14, 2026 and includes GPT-5-family forward-compatibility notes.

### MCP and A2A

- MCP intro: https://modelcontextprotocol.io/introduction
- MCP spec: https://modelcontextprotocol.io/specification/2025-11-25/basic/index
- MCP servers repo: https://github.com/modelcontextprotocol/servers
- MCP registry repo: https://github.com/modelcontextprotocol/registry
- A2A docs: https://a2aproject.github.io/A2A/latest/
- A2A Python SDK: https://github.com/a2aproject/a2a-python
- A2A JavaScript SDK: https://github.com/a2aproject/a2a-js

## Resulting repository decisions

1. Keep OpenAI Swarm only as transition material.
2. Add a dedicated agent evals and observability skill.
3. Prefer official releases pages in source maps for fast-moving frameworks.
4. Keep browser-native and local-first agent surfaces clearly separated.
5. Keep model-pinned frontier examples only as migration references, not as first-choice operational skills.
