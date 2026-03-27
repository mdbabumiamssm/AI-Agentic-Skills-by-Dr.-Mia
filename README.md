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

- OpenAI
- Anthropic
- Google Gemini
- Amazon Bedrock
- Azure AI Foundry
- Cohere
- Mistral
- DeepSeek
- xAI Grok
- provider and framework GitHub maintenance audits

### 2. Platform and orchestration code
The repository is not documentation-only. It also includes platform code for kernels, adapters, evaluators, reports, schemas, and utilities under `platform/`, plus supporting code under `src/`.

### 3. Research and analytics assets
The repository includes reference material and tutorial assets for practical learning and experimentation, especially under `Skills/Research_Tools/Data_Analysis/tutorials/`, covering:

- Python analytics workflows
- R tidyverse training material
- SQL datasets and learning assets
- Tableau and Power BI visualization content

These tutorial directories are now tracked directly inside the parent repository, which makes the repository easier to version, review, and distribute.

## Best entry points

If you are new to the repository, start here:

| Goal | Where to start |
|---|---|
| Explore domain skills | `Skills/` |
| Work on orchestration or kernel logic | `platform/` |
| Inspect supporting agent code | `src/` |
| Review tests and sample outputs | `tests/`, `test_demonstration/` |
| Browse research and planning material | `docs/` |
| Use tutorial/reference datasets | `Skills/Research_Tools/Data_Analysis/tutorials/` |

## Repository structure

```text
AIAGENTICSKILLS/
|-- Skills/                # Primary skill library across AI, medicine, science, law, finance, and engineering
|-- platform/              # Kernel, adapters, evaluator, reports, schema, examples, utilities
|-- src/                   # Supporting source code for research and chemistry agents
|-- docs/                  # Research notes, strategy docs, presentations
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

Automation:

- use `bash scripts/run_codex_repo_refresh.sh` for a Codex-driven provider/framework refresh workflow
- see `automation/CODEX_REPO_REFRESH.md` for run modes, prerequisites, and publication options

## Design principles

This repository is organized around practical reuse.

- skills are grouped by domain rather than by model vendor alone
- research assets and implementation assets coexist to reduce context switching
- platform code is separated from skill content to keep orchestration maintainable
- tutorial collections are kept inside the repository when they are needed as first-class reference material

## Contributing

Contributions should keep the repository navigable and self-contained.

Recommended practice:

1. create a focused feature branch
2. update the relevant skill, platform, or documentation area
3. keep `README.md` and related docs aligned with structural changes
4. run any relevant validation for the component you changed
5. open a pull request with a concise technical summary

## License

Copyright (c) 2026 MD BABU MIA, PhD.
All rights reserved.

This repository includes open-source components and reference material under their respective licenses. Repository-specific curation, structure, and workflow organization remain attributable to the author.

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
