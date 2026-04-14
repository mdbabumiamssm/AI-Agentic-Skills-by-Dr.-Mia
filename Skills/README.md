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

# Skills Repository (2026 Edition)

> **A production-focused repository for biomedical, clinical, and agentic AI systems.**

![Status](https://img.shields.io/badge/Status-Active-green)
![Agents](https://img.shields.io/badge/Agents-Orchestrated-blue)
![Domain](https://img.shields.io/badge/Domain-Biotech%20%7C%20Clinical%20%7C%20Genomics-purple)
![Tech](https://img.shields.io/badge/Tech-MCP%20%7C%20A2A%20%7C%20Frameworks-orange)

## Overview

This repository is a comprehensive library of skills, agents, framework operations, mathematical foundations, and domain workflows for modern (2026) AI systems. Unlike standard chatbot repos, this project is organized around agentic workflows that plan, execute, use tools, maintain state, and coordinate across specialized components.

The current curation standard is practical rather than hype-driven:

* use official docs and official GitHub repositories for operational facts
* favor maintained frameworks over abandoned demos
* separate protocol interoperability from provider-specific abstractions
* keep domain workflows adjacent to the framework guidance needed to run them

## Current Agentic AI Focus (checked April 7, 2026)

This collection now emphasizes the stacks that are most relevant for production-oriented agent development:

* **OpenAI Agent Platform:** Responses API, Agents SDK, AgentKit, built-in tools, evals, and deployment patterns.
* **LangGraph + LangSmith:** Durable execution, stateful graphs, interrupts, and observability for long-running agents.
* **PydanticAI + Logfire:** Typed Python agents, schema-first outputs, provider portability, and operational tracing.
* **Google ADK + A2A:** Code-first Google agent development with explicit multi-agent interoperability paths.
* **Microsoft Agent Framework:** Unified Python/.NET path for teams consolidating prior AutoGen and Semantic Kernel estates.
* **Agno:** Unified Python runtime for agents, teams, workflows, knowledge, and production serving.
* **CrewAI:** Role-based agent crews plus event-driven flows with first-party tracing.
* **MCP + A2A:** Clear separation between tool/resource interoperability and agent-to-agent interoperability.

## Key Capabilities

### Agentic AI

* **Framework Selection:** `Agentic_AI/AI_Framework_GitHub_Maintainers_2026/` tracks official docs and repos for maintained framework choices.
* **LangGraph Operations:** `Agentic_AI/LangGraph_Production_2026/` focuses on checkpoints, interrupts, human review, and durable orchestration.
* **PydanticAI Operations:** `Agentic_AI/PydanticAI_2026/` focuses on typed agents, structured outputs, tool schemas, and UI/event streaming.
* **Google ADK:** `Agentic_AI/Google_ADK_2026/` captures ADK quickstarts, evaluation, deployment, and A2A-aware design.
* **Microsoft Agent Framework:** `Agentic_AI/Microsoft_Agent_Framework_2026/` covers graph workflows, memory, and migration direction.
* **Agno:** `Agentic_AI/Agno_Operations_2026/` covers agents, teams, workflows, knowledge, approvals, and AgentOS deployment.
* **CrewAI:** `Agentic_AI/CrewAI_Production_2026/` covers crews, flows, memory, and tracing.
* **A2A Protocol:** `Agentic_AI/A2A_Protocol_2026/` documents how independently hosted agents should interoperate.
* **MCP Tool Surfaces:** `Agentic_AI/MCP_Registry/` documents tool discovery and shared server patterns.
* **Architecture Patterns:** `Agentic_AI/Agent_Architectures/`, `Agentic_AI/Multi_Agent_Systems/`, and `Agentic_AI/Memory_Systems/` keep reusable coordination patterns close to the framework layer.

### AI Providers

Current first-party AI platform coverage includes dedicated 2026 operational skills for:

* OpenAI
* Anthropic
* Google Gemini
* Amazon Bedrock
* Azure AI Foundry
* Cohere
* Mistral
* DeepSeek
* xAI Grok
* provider and framework GitHub maintenance audits

### Domain Skill Libraries

The broader `Skills/` tree continues to cover:

* clinical and healthcare automation
* genomics, transcriptomics, epigenomics, and multi-omics
* drug discovery and structural biology
* lab automation and self-driving labs
* software engineering and analytics
* finance, legal, and cross-domain productivity workflows

## Directory Structure

```text
Skills/
├── Agentic_AI/           # Framework operations, protocols, orchestration patterns, memory, web agents
├── AI_Providers/         # First-party provider operations and provider maintenance guidance
├── Clinical/             # Clinical workflows, trial matching, prior auth, pathology, radiology
├── Drug_Discovery/       # Molecular design, structure, tox, chemistry tools
├── Genomics/             # Single cell, spatial, CRISPR, variant analysis, long-read workflows
├── Research_Tools/       # Literature, observability, reporting, search, pathway analysis
├── Software_Engineering/ # Core Python, web development, codebase investigation
└── ...                   # Additional domain libraries across science, medicine, finance, and engineering
```

## How to Start

If your entry point is agentic AI rather than a scientific domain:

1. read `Agentic_AI/AI_Framework_GitHub_Maintainers_2026/`
2. choose the framework skill that matches your constraints
3. decide whether you need MCP, A2A, or both
4. only then move into the domain skill directories

For a curated map of the current official docs, repos, and benchmark papers that should anchor future refreshes, read `docs/strategy/AGENTIC_AI_RESOURCE_MAP_2026_04.md`.

This prevents starting from the wrong orchestration stack and then forcing domain workflows to fit it later.

## Roadmap Direction

The repository's current direction is to strengthen:

* durable framework guidance
* protocol interoperability
* operational observability and evals
* framework selection logic for domain teams
* reusable skill formats that stay current as provider ecosystems evolve

The repository should no longer present older agent patterns as if they are sufficient by themselves for production use.

---
*Maintained by the Artificial Intelligence Group.*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
