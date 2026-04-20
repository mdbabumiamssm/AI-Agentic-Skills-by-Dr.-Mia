# AI Providers (2026)

This directory is the repository's provider and delivery-surface layer: first-party API operations, provider SDK migration guidance, coding-agent surfaces, ChatGPT app-builder patterns, and maintenance health signals grounded in official docs and official repositories.

## Current Priority Areas (checked April 20, 2026)

### 1. Core Provider Operations

Use these when you are choosing models, SDKs, auth paths, or deployment boundaries for application workloads:

- `OpenAI_Platform_Operations_2026/`
- `OpenAI_Agents_SDK_2026/`
- `Anthropic_Claude_Operations_2026/`
- `Google_Gemini_Operations_2026/`
- `AWS_Bedrock_Operations_2026/`
- `Azure_AI_Foundry_Operations_2026/`
- `Cloud_AI_Operations_AWS_Azure_2026/`
- `Cohere_Platform_Operations_2026/`
- `Mistral_Platform_Operations_2026/`
- `DeepSeek_API_Operations_2026/`
- `xAI_Grok_Operations_2026/`
- `Frontier_OSS_Models_2026/`

### 2. Coding-Agent and Developer Surfaces

Use these when the primary requirement is code execution, repository automation, developer workflow acceleration, or ChatGPT-native app delivery rather than a generic runtime abstraction:

- `OpenAI_Codex_Operations_2026/` for Codex app, CLI, IDE, SDK, AGENTS.md, MCP, skills, plugins, subagents, and GitHub automation.
- `Claude_Code_Operations_2026/` for Anthropic's terminal-native coding agent, agent teams, Agent SDK, and GitHub Actions workflows.
- `Gemini_CLI_Operations_2026/` for Google's terminal-native Gemini agent, MCP integration, grounding, checkpointing, and headless automation.
- `OpenAI_Apps_SDK_2026/` for building ChatGPT-native MCP apps, widgets, metadata, and submission-ready integrations.

These are not replacements for framework selection. They are delivery surfaces and operator environments.

### 3. Maintenance and Selection Guidance

- `AI_Provider_GitHub_Maintainers_2026/` for official repository freshness, release cadence, and current maintenance signals.

Use this before adopting community wrappers or old examples that have drifted away from the official SDK and docs path.

## How to Choose

### Choose a core provider skill when:

- you are wiring models, auth, SDKs, or deployment targets into an application
- you need migration guidance for changing model families or deprecated endpoints
- you are deciding between direct provider APIs and managed clouds such as Bedrock or Azure AI Foundry

### Choose a coding-agent or app-builder skill when:

- the main deliverable is code change automation, CI-driven repository work, or terminal-native coding assistance
- you need repo-local instruction files such as `AGENTS.md`, `CLAUDE.md`, or `GEMINI.md`
- you need ChatGPT app distribution via MCP-backed tools and widgets

### Pair this directory with `Skills/Agentic_AI/` when:

- you also need durable orchestration, graph state, typed multi-agent workflows, or explicit eval pipelines
- the coding agent is only one execution surface inside a broader agent system

## Repository Curation Standard

Provider content in this directory should meet five bars:

1. Official documentation and official GitHub repositories are the primary operational sources.
2. Delivery surfaces such as Codex, Claude Code, Gemini CLI, and Apps SDK stay separate from generic framework guidance.
3. Model recommendations, auth assumptions, and deployment boundaries are stated explicitly.
4. Cloud-provider wrappers are documented as distinct from direct provider SDK paths.
5. Repo-maintenance health should come from official release and repository activity, not secondary blog summaries.

## Recommended Reading Order

1. `AI_Provider_GitHub_Maintainers_2026/`
2. `OpenAI_Platform_Operations_2026/`
3. `Anthropic_Claude_Operations_2026/`
4. `Google_Gemini_Operations_2026/`
5. `AWS_Bedrock_Operations_2026/`
6. `Azure_AI_Foundry_Operations_2026/`
7. `Cloud_AI_Operations_AWS_Azure_2026/`
8. `OpenAI_Codex_Operations_2026/`
9. `Claude_Code_Operations_2026/`
10. `Gemini_CLI_Operations_2026/`
11. `OpenAI_Apps_SDK_2026/`
12. `Frontier_OSS_Models_2026/`

## Current Direction

The repository should keep strengthening:

- first-party API operations anchored in current docs
- coding-agent workflows that reflect the real tools developers now use
- ChatGPT-native app and MCP delivery guidance
- security-aware local, cloud, and CI execution boundaries
- selection logic that keeps provider, framework, and protocol decisions clearly separated
