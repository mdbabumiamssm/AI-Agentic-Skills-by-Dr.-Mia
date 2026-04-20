# AI-Agentic-Skills-by-Dr.-Mia

<p align="left">
  <img src="https://img.shields.io/badge/Status-Active-green" alt="status" />
  <img src="https://img.shields.io/badge/Collection-AI%20Agentic%20Skills-blue" alt="collection" />
  <img src="https://img.shields.io/badge/Focus-Biomed%20%7C%20AI%20%7C%20Coding%20%7C%20Analytics-orange" alt="focus" />
  <img src="https://img.shields.io/badge/License-MIT-yellow" alt="license" />
</p>

A curated, multi-domain repository of agentic AI skills, workflow definitions, platform components, research utilities, and applied tutorial assets for building real-world AI systems in biomedical science, healthcare, software engineering, analytics, automation, and decision support.

⚠️ IMPORTANT DISCLAIMER & COPYRIGHT NOTICE

This repository, its architecture, agent designs, and specific implementations are the intellectual property of MD BABU MIA, PhD.

While open-source components are licensed under MIT, the unique curation, architecture, and agentic workflows are proprietary to the author.

If you fork, clone, or copy this repository for public use, you MUST:

Retain this copyright notice.
Explicitly credit MD BABU MIA, PhD as the original author.
Link back to the original repository.
Plagiarism or uncredited redistribution is strictly prohibited.

## Why this repository exists

Most AI repositories are either model demos, narrow toolkits, or disconnected prompt libraries. This repository is structured differently: it is a working skill system.

It brings together:

- reusable skill definitions under `Skills/`
- orchestration and kernel components under `platform/`
- supporting agent code under `src/`
- validation material under `tests/` and `test_demonstration/`
- documentation and strategy material under `docs/`
- directly tracked analytics tutorial collections under `Skills/Research_Tools/Data_Analysis/tutorials/`
- repository-local research snapshots under `sources/` when a refresh depends on current upstream verification

The result is a practical foundation for teams building AI agents that need domain knowledge, repeatable workflows, implementation scaffolding, and real reference assets in one place.

## What you will find here

### 1. Domain-specific skill libraries
The `Skills/` tree contains domain-organized assets for agent behavior, tool usage, research workflows, analysis playbooks, and implementation references.

Representative areas include:

- `Skills/Agentic_AI`
- `Skills/AI_Providers`
- `Skills/Software_Engineering`
- `Skills/Research_Tools`
- `Skills/Clinical`
- `Skills/Genomics`
- `Skills/Drug_Discovery`
- `Skills/Finance`
- `Skills/Legal`
- `Skills/Science`
- `Skills/User_Collections`

Additional specialized coverage spans pathology AI, radiology AI, oncology, precision medicine, protein science, structural biology, synthetic biology, metabolomics, microbiome, workflow management, and more.

Current first-party AI platform coverage includes dedicated 2026 operational skills for:

- OpenAI platform operations
- OpenAI Agents SDK
- OpenAI Codex operations
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

Current agent-framework and interoperability coverage now includes dedicated 2026 skills for:

- LangGraph production patterns
- PydanticAI
- Google Agent Development Kit (ADK)
- Microsoft Agent Framework
- Agno operations
- CrewAI production patterns
- Mastra production patterns
- LlamaIndex workflows
- smolagents
- DSPy optimization workflows
- Agent2Agent (A2A) protocol
- Model Context Protocol (MCP) registry and server patterns
- Browser Use for browser-native execution
- agent evals and observability

### 2. Platform and orchestration code
The repository is not documentation-only. It also includes platform code for kernels, adapters, evaluators, reports, schemas, and utilities under `platform/`, plus supporting code under `src/`.

### 3. Research and analytics assets
The repository includes reference material and tutorial assets for practical learning and experimentation, especially under `Skills/Research_Tools/Data_Analysis/tutorials/`, covering:

- Python analytics workflows
- R tidyverse training material
- SQL datasets and learning assets
- Tableau and Power BI visualization content

These tutorial directories are tracked directly inside the parent repository, which makes the repository easier to version, review, and distribute.

## Best entry points

If you are new to the repository, start here:

| Goal | Where to start |
|---|---|
| Explore domain skills | `Skills/` |
| Choose a production agent stack | `Skills/Agentic_AI/AI_Framework_GitHub_Maintainers_2026/` |
| Choose a provider runtime | `Skills/AI_Providers/` |
| Review coding-agent and app-builder surfaces | `Skills/AI_Providers/README.md` |
| Review current official source maps | `docs/strategy/AGENTIC_AI_RESOURCE_MAP_2026_04.md` |
| Review upstream maintenance notes | `docs/strategy/AGENTIC_AI_UPSTREAM_HEALTH_2026_04.md` |
| Inspect supporting code | `platform/`, `src/` |
| Review tests and sample outputs | `tests/`, `test_demonstration/` |
| Inspect saved refresh research | `sources/` |

## Repository structure

```text
AIAGENTICSKILLS/
|-- Skills/                # Primary skill library across AI, medicine, science, law, finance, and engineering
|-- platform/              # Kernel, adapters, evaluator, reports, schema, examples, utilities
|-- src/                   # Supporting source code for research and chemistry agents
|-- docs/                  # Research notes, strategy docs, presentations
|-- sources/               # Saved current-source research snapshots used during refreshes
|-- tests/                 # Validation assets and test material
|-- test_demonstration/    # Demonstration workflows and outputs
|-- archive/               # Archived material and prior snapshots
|-- figures/               # Visual assets and figures
|-- cache/                 # Cached artifacts
`-- README.md
```

## Typical use cases

This repository is useful if you are building or curating AI systems for:

- biomedical and translational research workflows
- clinical reasoning and healthcare-oriented automation
- scientific analysis and lab-support tooling
- software engineering agents and code workflows
- data analysis, reporting, and dashboard generation
- agentic orchestration across tools, providers, and execution environments

## Current curation priorities

As of April 20, 2026, the highest-priority repository upgrades are:

- first-party provider operations grounded in official docs
- maintained agent frameworks rather than stale demo stacks
- coding-agent and ChatGPT app-builder surfaces such as Codex, Claude Code, Gemini CLI, and Apps SDK
- explicit transition labeling for legacy patterns such as OpenAI Swarm and model-pinned frontier examples
- browser-native and local-first agent execution documented as distinct choices
- protocol interoperability through MCP and A2A
- evaluation, tracing, and observability as first-class design constraints
- source maps and research snapshots that make refreshes reproducible
- strategy documents that explain which framework to choose and why

## Quick start

Clone the repository:

```bash
git clone https://github.com/mdbabumiamssm/AI-Agentic-Skills-by-Dr.-Mia.git
cd AI-Agentic-Skills-by-Dr.-Mia
```

This repository contains multiple independent areas rather than one single runtime. Install dependencies according to the part of the repository you are using.

Typical workflow:

1. identify the domain area you need under `Skills/`
2. inspect the related `SKILL.md` and implementation assets
3. use `platform/` if you need orchestration, evaluation, or kernel support
4. validate with the relevant test or demonstration material
5. consult `docs/strategy/` and `sources/` when the work depends on current provider or framework state

Automation:

- use `bash scripts/run_codex_repo_refresh.sh` for a Codex-driven provider and framework refresh workflow
- see `automation/CODEX_REPO_REFRESH.md` for run modes, prerequisites, and publication options

## Design principles

This repository is organized around practical reuse.

- skills are grouped by domain rather than by model vendor alone
- official docs and official repositories are the primary operational anchors
- transition material is kept, but it is labeled clearly so it does not masquerade as a default recommendation
- platform code is separated from skill content to keep orchestration maintainable
- repository-local research snapshots are saved when current upstream verification affects curation decisions

## Contributing

Contributions should keep the repository navigable and current.

Recommended practice:

1. create a focused feature branch
2. update the relevant skill, platform, or documentation area
3. keep `README.md`, `Skills/README.md`, and strategy docs aligned with structural changes
4. run the repository refresh validation on changed or affected skills
5. open a pull request with a concise technical summary

## License

Copyright (c) 2026 MD BABU MIA, PhD.
All rights reserved.

This repository includes open-source components and reference material under their respective licenses. Repository-specific curation, structure, and workflow organization remain attributable to the author.

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
