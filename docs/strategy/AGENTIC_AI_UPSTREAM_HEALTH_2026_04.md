# Agentic AI Upstream Health (April 16, 2026)

This note captures the current upstream health signals that matter most for repository curation. The goal is not to mirror every project on GitHub. The goal is to keep this repository anchored to the maintained stacks, the live interoperability standards, and the browser or local execution surfaces that still move the field forward.

## First-choice framework signals

### OpenAI Agents SDK

- Official docs remain the primary anchor for runtime selection, built-in tools, tracing, evals, and UI surfaces.
- The official Python repository is active and the releases page shows `v0.14.1` published on April 15, 2026.
- Current guidance should prefer the Agents SDK over Swarm for production work.

### LangGraph

- LangGraph remains active, but versioning is split across packages. The releases page shows `langgraph-checkpoint==4.0.2` on April 15, 2026 and alpha releases for the main `langgraph` package on April 14, 2026.
- The checkpoint package is operationally important enough to track separately because durable execution is the core reason to choose LangGraph.

### PydanticAI

- The project is highly active. The releases page shows `v1.83.0` on April 16, 2026.
- Recent changes include FastMCP toolset capabilities and prompt caching integrations, which reinforces the value of keeping typed-tool and observability guidance current.

### Google ADK

- The repository is moving quickly. The stable line shows `v1.30.0` on April 13, 2026, while the same release page also exposes active `v2.0.0a3` pre-release work from April 9, 2026.
- The release notes show A2A, MCP toolset, live mode, environment tools, and skill-system work all changing at once, so model and workflow guidance should stay tightly pinned to official docs.

### Microsoft Agent Framework

- The release page shows the Python beta line at `1.0.0b260107` on January 7, 2026.
- It remains important for teams consolidating AutoGen and Semantic Kernel estates, but the beta status means repository guidance should keep rollout and compatibility notes close to the sources.

## Strong-secondary framework signals

### CrewAI

- The official releases page shows `1.8.0` on January 8, 2026.
- Recent features include A2A hooks, human-in-the-loop flow controls, and stronger event handling, which keeps CrewAI relevant for deliberate role-based designs.

### Agno

- The official releases page shows `v2.3.24` on January 8, 2026.
- The project remains useful as an integrated runtime spanning agents, workflows, tools, and knowledge, but should still be chosen for fit rather than as a default replacement for everything else.

## Browser and local execution signals

### Browser Use

- The official releases page shows `0.12.6` on April 2, 2026.
- The same release stream shows a March 25, 2026 security fix removing `litellm` from core dependencies after a supply-chain incident, which is a reminder to track browser-agent dependency churn directly from upstream.
- Browser Use remains a strong fit when managed browser infrastructure and persistent browser sessions are central.

### OpenClaw

- The official releases page shows `openclaw 2026.4.14` on April 14, 2026.
- Current release notes explicitly mention forward-compatibility work for GPT-5 family models, which is a useful signal that OpenClaw remains active as a local automation surface.
- OpenClaw should be treated as local-first automation, not as a replacement for managed browser infrastructure.

## Coding-agent and app-builder surfaces

### OpenAI Codex

- The official Codex docs now expose a broad operational surface: app, IDE extension, CLI, web, integrations, AGENTS.md, MCP, plugins, skills, subagents, SDK, app server, MCP server, and GitHub Action.
- The current OpenAI models catalog recommends GPT-5.4 as the default starting point for complex reasoning and coding, with GPT-5.4 mini explicitly positioned for coding, computer use, and subagents.
- Repository guidance should treat Codex as a first-class coding-agent surface, not as a hidden footnote under general OpenAI API operations.

### Claude Code

- The official overview positions Claude Code as Anthropic's terminal-native coding tool and now highlights multiple Claude Code agents, the Agent SDK, GitHub Actions, GitLab CI/CD, channels, routines, and Chrome debugging.
- The official `anthropics/claude-code` repository currently deprecates npm installation in favor of install scripts, Homebrew, and Windows-native installers, which is operationally important enough to capture explicitly in skill guidance.
- The official `anthropics/claude-code-action` repository remains the primary automation surface for GitHub-native repository workflows.

### Gemini CLI

- The official `google-gemini/gemini-cli` repository is highly active and currently shows `v0.38.2` released on April 17, 2026.
- The repository positions Gemini CLI as an open-source terminal AI agent with MCP integration, Google Search grounding, checkpointing, headless scripting, GitHub integration, and `GEMINI.md` context files.
- Weekly stable and preview release channels plus nightly builds make release cadence itself a relevant operational signal for curation.

### OpenAI Apps SDK

- The official Apps SDK docs now define a full ChatGPT app lifecycle around MCP apps, server setup, UI widgets, metadata optimization, deployment, testing, submission, and reference material.
- The Apps SDK home page still states that the current public distribution path is submission-based, while approved apps appear in the ChatGPT apps store and OpenAI creates a plugin for Codex distribution.
- The official examples repository remains the best starting point for combining MCP servers, structured tool results, and widget UI metadata.

## Interoperability standards

### MCP

- The Model Context Protocol remains the default standard for exposing tools and resources to agents.
- Repository guidance should keep using official MCP documentation and official `modelcontextprotocol` repositories as the source of truth.

### A2A

- Agent2Agent remains the right standard to track when independently hosted agents must communicate across framework boundaries.
- It should remain clearly separated from MCP in repository docs and skills.

## Transition material

### OpenAI Swarm

- The official `openai/swarm` repository now states that Swarm is experimental and educational, has been replaced by the OpenAI Agents SDK, and should not be used for new production builds.
- It still has value as a handoff teaching reference, but not as a first-choice stack.

## Curation actions implied by these signals

1. Keep first-choice framework skills near official docs and official repositories.
2. Preserve legacy or model-pinned examples only when they are clearly marked as transition material.
3. Add observability and evaluation guidance as a first-class repository skill instead of burying it inside framework notes.
4. Track security-sensitive browser and local automation stacks through official release pages, not blog summaries.
5. Keep coding-agent and ChatGPT app-builder surfaces as first-class provider skills rather than burying them in generic framework documentation.
6. Keep resource maps opinionated: fewer stronger anchors are better than many overlapping references.
