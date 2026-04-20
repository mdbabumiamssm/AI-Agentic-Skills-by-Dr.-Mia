# Agentic AI (2026)

This directory is the repository's control plane for modern agent systems: orchestration patterns, framework operations, protocol interoperability, memory, browser execution, local automation, and multi-agent coordination.

## Current Priority Areas (checked April 16, 2026)

### 1. Production Agent Frameworks

Use these directories first when choosing a framework for new builds:

- `LangGraph_Production_2026/` for durable execution, checkpoints, interrupts, and long-running stateful workflows.
- `PydanticAI_2026/` for typed Python agents, structured outputs, and schema-first service boundaries.
- `Google_ADK_2026/` for Google ecosystem deployments, evaluation, and A2A-aware orchestration.
- `Microsoft_Agent_Framework_2026/` for Python or .NET enterprise workflows and migrations away from split AutoGen and Semantic Kernel estates.

### 2. Strong Secondary Frameworks

Use these when their fit is deliberate rather than assumed:

- `Agno_Operations_2026/` for one Python runtime spanning agents, teams, workflows, knowledge, approvals, and served production APIs.
- `CrewAI_Production_2026/` for role-based crews plus explicit event-driven flows.
- `Mastra_Production_2026/` for TypeScript-first product teams that want agents, workflows, MCP, evals, and observability together.
- `LlamaIndex_Workflows_2026/` for document- and data-aware agent workflows with explicit events and serving paths.
- `Smolagents_2026/` for lightweight Python code agents and tool-calling agents with minimal abstraction overhead.

### 3. Specialized Program Layers

- `DSPy_2026/` for modular LM programming, prompt optimization, compilation, and evaluation-first improvement loops.
- `Agent_Evals_Observability_2026/` for offline datasets, rollout gates, traced runs, and production monitoring.

These are not standalone reasons to choose a whole runtime. They are layers that help make the chosen runtime reliable.

### 4. Coding-Agent and App Surfaces

The repository now treats coding agents and ChatGPT app-builder surfaces as a separate provider-layer concern rather than burying them inside framework docs:

- see `../AI_Providers/OpenAI_Codex_Operations_2026/`
- see `../AI_Providers/Claude_Code_Operations_2026/`
- see `../AI_Providers/Gemini_CLI_Operations_2026/`
- see `../AI_Providers/OpenAI_Apps_SDK_2026/`

These are execution and delivery surfaces. They are useful when the main job is code automation, CI integration, or ChatGPT-native UI delivery. They should complement, not replace, deliberate framework selection.

### 5. Protocol Interoperability

- `A2A_Protocol_2026/` for agent-to-agent interoperability across framework boundaries.
- `MCP_Registry/` for Model Context Protocol discovery, tool exposure, and shared tool surfaces.

These solve different problems:

- **MCP** standardizes agent-to-tool and agent-to-resource access.
- **A2A** standardizes communication between remote agents as networked services.

### 6. Browser and Local Execution

- `Web_Agents/Browser_Use_2026/` for browser-native execution when live web work, persistent browser workspaces, and managed infrastructure are core to the system.
- `Web_Agents/OpenClaw_Agent/` for local-first automation when browser work is only one part of a broader machine-local execution surface.

Browser execution should usually be treated as a substrate under a larger orchestrator, not as a replacement for framework selection.

### 7. Reusable Architecture Patterns

- `Multi_Agent_Systems/` for supervisor-worker orchestration and debate patterns.
- `Agent_Architectures/` for ReAct, plan-and-solve, and self-correction references.
- `Memory_Systems/` for short-term, episodic, and persistent memory patterns.
- `Frontier_Models/` for transition references, not default operational guidance.

### 8. Transition and Legacy Material

- `OpenAI_Swarm/` remains useful only as a lightweight historical reference and migration bridge.
- `Frontier_Models/` contains model-pinned examples that should be translated into current provider skills before real deployment.
- Framework selection and maintenance decisions should be grounded in `AI_Framework_GitHub_Maintainers_2026/`.

For greenfield builds, do not start from older assumptions that every agent system should be a prompt loop plus ad hoc tool calls. Prefer framework-managed state, typed interfaces, explicit evals, and interoperable protocols.

## Repository Curation Standard

Agentic AI content in this repository should now meet five bars:

1. Official documentation and official GitHub repositories are the primary sources for operational facts.
2. New framework skills must explain where they are strongest and where they should not be the default.
3. Interoperability guidance must distinguish between MCP, A2A, provider-native tools, and framework-specific abstractions.
4. Index files should surface current first-choice stacks instead of leaving them buried in subdirectories.
5. Legacy or model-pinned material must be labeled as transition guidance rather than current best practice.

## Recommended Reading Order

1. `AI_Framework_GitHub_Maintainers_2026/`
2. `LangGraph_Production_2026/`
3. `PydanticAI_2026/`
4. `Google_ADK_2026/`
5. `Microsoft_Agent_Framework_2026/`
6. `Agent_Evals_Observability_2026/`
7. `Mastra_Production_2026/`
8. `LlamaIndex_Workflows_2026/`
9. `Smolagents_2026/`
10. `DSPy_2026/`
11. `Agno_Operations_2026/`
12. `CrewAI_Production_2026/`
13. `A2A_Protocol_2026/`
14. `MCP_Registry/`
15. `Web_Agents/Browser_Use_2026/`
16. `Web_Agents/OpenClaw_Agent/`
17. `OpenAI_Swarm/` only if you are translating legacy material

## Previous and Ongoing Work

The existing code modules in this directory remain useful:

- `Multi_Agent_Systems/orchestrator.py`
- `Multi_Agent_Systems/debate_supervisor.py`
- `Agent_Architectures/Plan_and_Solve/`
- `Agent_Architectures/ReAct_Agent/`
- `Memory_Systems/memory_architecture.py`

These are pattern references. The 2026 framework, protocol, evaluation, and execution skills are the operational layer that should guide how new systems are actually built.
