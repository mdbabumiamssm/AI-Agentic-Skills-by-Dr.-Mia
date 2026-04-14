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

# Agentic AI Stack Blueprint (April 13, 2026)

This document defines the recommended agentic AI stack priorities for this repository as of April 7, 2026.

## What this repository should optimize for

The repository should prefer:

1. maintained first-party or strongly maintained framework ecosystems
2. protocol interoperability instead of one-off integration layers
3. typed interfaces, evals, and observability instead of prompt-only orchestration
4. domain workflows that can be ported across providers when needed

## Recommended framework choices by use case

### Use `OpenAI_Agents_SDK_2026` when:

- the workload is already standardized on OpenAI
- built-in hosted tools matter
- you want a fast path to tracing, evals, UI embeds, or background execution

### Use `LangGraph_Production_2026` when:

- you need durable execution and resumability
- long-running multi-step workflows are central
- human review and interrupt/resume control are required
- state transitions must be explicit

### Use `PydanticAI_2026` when:

- typed Python interfaces and structured outputs are critical
- the boundary between agent calls and application code should stay explicit
- you want clean schema enforcement and provider portability

### Use `Google_ADK_2026` when:

- the deployment target leans into Google tooling
- you want code-first orchestration with evaluation support
- A2A interoperability is a design requirement

### Use `Mastra_Production_2026` when:

- the team is TypeScript-first
- product-facing agents, workflows, MCP, evals, and observability should live in one stack
- the agent runtime is part of the application, not a separate research sandbox

### Use `LlamaIndex_Workflows_2026` when:

- the core problem is document, data, or retrieval-aware orchestration
- explicit events and workflow serving matter
- the agent should stay close to indexing and retrieval boundaries

### Use `Smolagents_2026` when:

- minimal Python abstractions are a feature, not a bug
- code agents or lightweight tool-calling agents fit better than a larger framework runtime
- the team wants flexible model backends with explicit sandbox choices

### Use `DSPy_2026` when:

- the hard problem is LM optimization, evaluation, or compilation
- prompt engineering needs to become modular and metric-driven
- the runtime framework can stay separate from the program optimization layer

### Use `Microsoft_Agent_Framework_2026` when:

- the team is split across Python and .NET
- there is existing AutoGen or Semantic Kernel history
- workflow graphs, threads, and enterprise hosting patterns matter

### Use `A2A_Protocol_2026` when:

- independently hosted agents need to talk to each other
- you want an inter-agent protocol rather than a framework-specific handoff API

### Use `MCP` when:

- agents need shared access to tools, data sources, or resources
- the same tool surfaces should be reusable across clients

## Frameworks to keep, but not treat as default greenfield choices

### Agno

Keep for:

- Python teams that want agents, teams, workflows, knowledge, and served runtime patterns in one stack
- systems that need approval workflows, session state, and production serving without stitching many layers

Do not treat as the first recommendation when a simpler first-tier choice already fits and the extra runtime surface would add unnecessary platform weight.

### CrewAI

Keep for:

- role-based multi-agent automation where crew and flow abstractions are a natural fit
- teams that want a Python-native framework with explicit multi-agent orchestration and first-party tracing

Do not treat as the universal default for all agent systems; choose it when crew/flow concepts are a real architectural benefit.

### Mastra

Keep for:

- TypeScript-first product teams that want agents, workflows, MCP, evals, and observability in one framework
- teams that need a direct path from local development to deployed product runtime

Do not treat as the first recommendation when the team is primarily Python-based or when a simpler provider-native approach is enough.

### LlamaIndex Workflows

Keep for:

- document, data, and retrieval-centric agent systems
- workflows that should stay close to indexing, ingestion, and serving concerns

Do not treat as the default orchestration layer for non-data-centric agent systems where another runtime is a better fit.

### smolagents

Keep for:

- lightweight code-agent or tool-calling-agent systems
- teams that want minimal abstraction overhead and explicit sandbox choices

Do not treat as the default when durable execution, human approval checkpoints, or large production control surfaces are central.

### DSPy

Keep for:

- modular LM programming, compilation, and evaluation-first optimization
- systems where metric-driven prompt and program improvement is central

Do not treat as the whole runtime stack when the real need is hosted orchestration, workflow state, or protocol interoperability.

### AutoGen

Keep for:

- teams already invested in AutoGen
- historical comparison and migration planning

Do not treat as the first recommendation for new Microsoft-centered builds without evaluating Microsoft Agent Framework first.

### Semantic Kernel

Keep for:

- existing production estates
- specific enterprise connectors or orchestration paths already built on it

Do not treat as the first greenfield default when Microsoft Agent Framework is a better fit.

### OpenAI Swarm

Keep for:

- small handoff demos
- historical reference

Do not treat as the primary production framework layer for new multi-agent systems.

## Protocol model

The repository should teach the following separation clearly:

- `MCP` = agent-to-tool and agent-to-resource interoperability
- `A2A` = agent-to-agent interoperability
- provider SDKs = provider-native inference, tools, and account-level capabilities
- frameworks = orchestration and application architecture

## High-value literature to keep close to the frameworks

Foundational agent papers:

- ReAct: https://arxiv.org/abs/2210.03629
- Toolformer: https://arxiv.org/abs/2302.04761
- Reflexion: https://arxiv.org/abs/2303.11366
- Tree of Thoughts: https://arxiv.org/abs/2305.10601
- Voyager: https://arxiv.org/abs/2305.16291
- Language Agent Tree Search: https://arxiv.org/abs/2310.04406
- GAIA benchmark: https://arxiv.org/abs/2311.12983
- SWE-agent: https://arxiv.org/abs/2405.15793
- DSPy: https://arxiv.org/abs/2310.03714
- OSWorld: https://os-world.github.io/

These are not substitutes for framework documentation, but they remain the right references for understanding why certain orchestration patterns exist.

## Official resources that should anchor future refreshes

- OpenAI Agents SDK: https://platform.openai.com/docs/guides/agents-sdk
- LangGraph: https://docs.langchain.com/oss/python/langgraph/overview
- LangSmith: https://docs.langchain.com/langsmith/observability
- PydanticAI: https://ai.pydantic.dev/
- Google ADK: https://google.github.io/adk-docs/
- Microsoft Agent Framework: https://learn.microsoft.com/en-us/agent-framework/get-started/
- Agno: https://docs.agno.com/
- CrewAI: https://docs.crewai.com/en
- Mastra: https://mastra.ai/docs
- LlamaIndex Workflows: https://developers.llamaindex.ai/python/llamaagents/workflows/
- smolagents: https://huggingface.co/docs/smolagents/en/index
- DSPy: https://dspy.ai/
- Browser Use: https://docs.browser-use.com/cloud/quickstart
- A2A Protocol: https://a2aproject.github.io/A2A/latest/
- Model Context Protocol: https://modelcontextprotocol.io/introduction

## Curation rules for future additions

- add a new framework skill only if it has active official docs and an official repository
- prefer a small number of high-signal skills over many shallow duplicates
- add selection logic, not just tool descriptions
- record official sources in `references/`
- update the index files when the recommendation set changes


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
