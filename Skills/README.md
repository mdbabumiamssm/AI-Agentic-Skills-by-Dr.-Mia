# Skills Repository (2026 Edition)

> **A production-focused repository for biomedical, clinical, and agentic AI systems.**

![Status](https://img.shields.io/badge/Status-Active-green)
![Agents](https://img.shields.io/badge/Agents-Orchestrated-blue)
![Domain](https://img.shields.io/badge/Domain-Biotech%20%7C%20Clinical%20%7C%20Genomics-purple)
![Tech](https://img.shields.io/badge/Tech-MCP%20%7C%20A2A%20%7C%20Frameworks-orange)

## Overview

This repository is a comprehensive library of skills, agents, framework operations, mathematical foundations, and domain workflows for modern (2026) AI systems. Unlike standard chatbot repos, this project is organized around agentic workflows that plan, execute, use tools, maintain state, and coordinate across specialized components.

The current curation standard is practical rather than hype-driven:

- use official docs and official GitHub repositories for operational facts
- favor maintained frameworks over abandoned demos
- separate protocol interoperability from provider-specific abstractions
- keep domain workflows adjacent to the framework guidance needed to run them
- preserve legacy material only when it is clearly marked as transition guidance

## Current Agentic AI Focus (checked April 20, 2026)

This collection emphasizes the stacks that are most relevant for production-oriented agent development:

- **OpenAI Agent Platform:** Responses API, Agents SDK, guardrails, built-in tools, tracing, evals, and deployment patterns.
- **LangGraph + LangSmith:** Durable execution, stateful graphs, interrupts, observability, and evaluation loops for long-running agents.
- **PydanticAI + Logfire:** Typed Python agents, schema-first outputs, provider portability, and operational tracing.
- **Google ADK + A2A:** Code-first Google agent development with explicit multi-agent interoperability paths.
- **Microsoft Agent Framework:** Unified Python and .NET path for teams consolidating prior AutoGen and Semantic Kernel estates.
- **Agno:** Unified Python runtime for agents, teams, workflows, knowledge, and production serving.
- **CrewAI:** Role-based agent crews plus event-driven flows with first-party observability.
- **MCP + A2A:** Clear separation between tool or resource interoperability and agent-to-agent interoperability.
- **Browser Use and OpenClaw:** Explicit separation between browser-native execution and local-first automation.
- **Agent Evals and Observability:** Regression loops, rollout gates, and live trace review as first-class design constraints.
- **Coding-Agent Surfaces:** OpenAI Codex, Claude Code, and Gemini CLI as first-class developer operating surfaces rather than side notes.
- **ChatGPT App Surfaces:** OpenAI Apps SDK for MCP-backed ChatGPT apps, widgets, metadata, and submission-ready integrations.

## Key Capabilities

### Agentic AI

- **Framework Selection:** `Agentic_AI/AI_Framework_GitHub_Maintainers_2026/` tracks official docs and repos for maintained framework choices.
- **LangGraph Operations:** `Agentic_AI/LangGraph_Production_2026/` focuses on checkpoints, interrupts, human review, and durable orchestration.
- **PydanticAI Operations:** `Agentic_AI/PydanticAI_2026/` focuses on typed agents, structured outputs, tool schemas, and UI or event streaming.
- **Google ADK:** `Agentic_AI/Google_ADK_2026/` captures ADK quickstarts, evaluation, deployment, and A2A-aware design.
- **Microsoft Agent Framework:** `Agentic_AI/Microsoft_Agent_Framework_2026/` covers graph workflows, memory, and migration direction.
- **Agno:** `Agentic_AI/Agno_Operations_2026/` covers agents, teams, workflows, knowledge, approvals, and AgentOS deployment.
- **CrewAI:** `Agentic_AI/CrewAI_Production_2026/` covers crews, flows, memory, and tracing.
- **A2A Protocol:** `Agentic_AI/A2A_Protocol_2026/` documents how independently hosted agents should interoperate.
- **MCP Tool Surfaces:** `Agentic_AI/MCP_Registry/` documents tool discovery and shared server patterns.
- **Browser and Local Execution:** `Agentic_AI/Web_Agents/` distinguishes Browser Use from OpenClaw instead of collapsing them into one category.
- **Evaluation and Tracing:** `Agentic_AI/Agent_Evals_Observability_2026/` centralizes rollout-safe agent quality practices.
- **Architecture Patterns:** `Agentic_AI/Agent_Architectures/`, `Agentic_AI/Multi_Agent_Systems/`, and `Agentic_AI/Memory_Systems/` keep reusable coordination patterns close to the framework layer.

### AI Providers

Current first-party AI platform coverage includes dedicated 2026 operational skills for:

- OpenAI platform operations
- OpenAI Agents SDK
- OpenAI Codex
- OpenAI Apps SDK
- Anthropic Claude
- Claude Code
- Google Gemini
- Gemini CLI
- Amazon Bedrock
- Azure AI Foundry
- Cohere
- Mistral
- DeepSeek
- xAI Grok
- provider and framework GitHub maintenance audits

### Domain Skill Libraries

The broader `Skills/` tree continues to cover:

- clinical and healthcare automation
- genomics, transcriptomics, epigenomics, and multi-omics
- drug discovery and structural biology
- lab automation and self-driving labs
- software engineering and analytics
- finance, legal, and cross-domain productivity workflows

## Directory Structure

```text
Skills/
├── Agentic_AI/           # Framework operations, protocols, orchestration patterns, memory, evals, browser and local execution
├── AI_Providers/         # First-party provider operations and provider maintenance guidance
├── Clinical/             # Clinical workflows, trial matching, prior auth, pathology, radiology
├── Drug_Discovery/       # Molecular design, structure, tox, chemistry tools
├── Genomics/             # Single cell, spatial, CRISPR, variant analysis, long-read workflows
├── Research_Tools/       # Literature, reporting, search, pathway analysis
├── Software_Engineering/ # Core Python, web development, codebase investigation
└── ...                   # Additional domain libraries across science, medicine, finance, and engineering
```

## How to Start

If your entry point is agentic AI rather than a scientific domain:

1. read `Agentic_AI/AI_Framework_GitHub_Maintainers_2026/`
2. review `AI_Providers/README.md` if your main need is a coding agent, ChatGPT app, or provider delivery surface
3. choose the framework skill that matches your constraints
4. decide whether you need MCP, A2A, browser-native execution, local-first automation, or some combination
5. add evaluation and observability requirements before broadening autonomy
6. only then move into the domain skill directories

For a curated map of the current official docs, repos, and benchmark papers that should anchor future refreshes, read `docs/strategy/AGENTIC_AI_RESOURCE_MAP_2026_04.md`.

For current upstream maintenance notes captured during the April 16 refresh, read `docs/strategy/AGENTIC_AI_UPSTREAM_HEALTH_2026_04.md` and `sources/research_20260416_agentic_ai_official_stack_health.md`.

## Roadmap Direction

The repository's current direction is to strengthen:

- durable framework guidance
- coding-agent and ChatGPT app-builder operations grounded in official docs
- protocol interoperability
- operational observability and evals
- framework selection logic for domain teams
- local and browser execution guidance with explicit security boundaries
- reusable skill formats that stay current as provider ecosystems evolve

The repository should no longer present older agent patterns as if they are sufficient by themselves for production use.
