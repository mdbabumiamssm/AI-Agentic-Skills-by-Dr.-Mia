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

# Agentic AI (2026)

This directory is the repository's control plane for modern agent systems: orchestration patterns, framework operations, protocol interoperability, memory, web execution, and multi-agent coordination.

## Current Priority Areas (checked April 13, 2026)

### 1. Production Agent Frameworks

Use these directories first when choosing a framework for new builds:

* `LangGraph_Production_2026/` for durable execution, checkpoints, interrupts, and long-running stateful workflows.
* `PydanticAI_2026/` for typed Python agents, structured outputs, and schema-first service boundaries.
* `Google_ADK_2026/` for Google ecosystem deployments, evaluation, and A2A-aware orchestration.
* `Microsoft_Agent_Framework_2026/` for Python/.NET enterprise workflows and migrations away from split AutoGen/Semantic Kernel estates.

### 2. Strong Secondary Frameworks

Use these when their fit is deliberate rather than assumed:

* `Agno_Operations_2026/` for one Python runtime spanning agents, teams, workflows, knowledge, approvals, and served production APIs.
* `CrewAI_Production_2026/` for role-based crews plus explicit event-driven flows.
* `Mastra_Production_2026/` for TypeScript-first product teams that want agents, workflows, MCP, evals, and observability together.
* `LlamaIndex_Workflows_2026/` for document and data-aware agent workflows with explicit events and serving paths.
* `Smolagents_2026/` for lightweight Python code agents and tool-calling agents with minimal abstraction overhead.

### 3. Specialized Program Layers

* `DSPy_2026/` for modular LM programming, prompt optimization, compilation, and evaluation-first improvement loops.

This is not the same thing as choosing the whole runtime. DSPy is strongest when the hard part is optimization and measurement rather than hosting.

### 4. Protocol Interoperability

* `A2A_Protocol_2026/` for agent-to-agent interoperability across framework boundaries.
* `MCP_Registry/` for Model Context Protocol discovery, tool exposure, and shared tool surfaces.

These solve different problems:

* **MCP** standardizes agent-to-tool and agent-to-resource access.
* **A2A** standardizes communication between remote agents as networked services.

### 5. Browser-Native Execution

* `Web_Agents/Browser_Use_2026/` for browser-native execution when live web work, persistent workspaces, and browser security boundaries are core to the system.

Browser execution should usually be treated as a substrate under a larger orchestrator, not as a replacement for framework selection.

### 6. Reusable Architecture Patterns

* `Multi_Agent_Systems/` for supervisor-worker orchestration and debate patterns.
* `Agent_Architectures/` for ReAct, plan-and-solve, and self-correction references.
* `Memory_Systems/` for short-term, episodic, and persistent memory patterns.
* `Web_Agents/` for browser-native execution paths.

### 7. Transition and Legacy Material

* `OpenAI_Swarm/` remains useful as a lightweight historical reference and for small handoff demos.
* Framework selection and maintenance decisions should be grounded in `AI_Framework_GitHub_Maintainers_2026/`.

For greenfield builds, do not start from older assumptions that every agent system should be a prompt loop plus ad hoc tool calls. Prefer framework-managed state, typed interfaces, explicit evals, and interoperable protocols.

## Repository Curation Standard

Agentic AI content in this repository should now meet four bars:

1. Official documentation and official GitHub repositories are the primary sources for operational facts.
2. New framework skills must explain where they are strongest and where they should not be the default.
3. Interoperability guidance must distinguish between MCP, A2A, provider-native tools, and framework-specific abstractions.
4. Index files should surface current first-choice stacks instead of leaving them buried in subdirectories.

## Recommended Reading Order

1. `AI_Framework_GitHub_Maintainers_2026/`
2. `LangGraph_Production_2026/`
3. `PydanticAI_2026/`
4. `Google_ADK_2026/`
5. `Microsoft_Agent_Framework_2026/`
6. `Mastra_Production_2026/`
7. `LlamaIndex_Workflows_2026/`
8. `Smolagents_2026/`
9. `DSPy_2026/`
10. `Agno_Operations_2026/`
11. `CrewAI_Production_2026/`
12. `A2A_Protocol_2026/`
13. `MCP_Registry/`
14. `Web_Agents/Browser_Use_2026/`

## Previous and Ongoing Work

The existing code modules in this directory remain useful:

* `Multi_Agent_Systems/orchestrator.py`
* `Multi_Agent_Systems/debate_supervisor.py`
* `Agent_Architectures/Plan_and_Solve/`
* `Agent_Architectures/ReAct_Agent/`
* `Memory_Systems/memory_architecture.py`

These are pattern references. The new 2026 framework and protocol skills are the operational layer that should guide how new systems are actually built.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
