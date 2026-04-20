# Coding-Agent and App Surface Research (April 20, 2026)

This note captures the official-source findings that justify adding first-class provider skills for coding agents and ChatGPT-native app delivery.

## OpenAI Codex

- The official Codex docs now expose a broad operator surface covering app, IDE extension, CLI, web, integrations, config, `AGENTS.md`, MCP, plugins, skills, subagents, SDK, app server, MCP server, and GitHub Action.
- The current OpenAI models catalog recommends GPT-5.4 as the default starting point for complex reasoning and coding, with GPT-5.4 mini explicitly positioned for coding, computer use, and subagents.
- This is distinct enough from generic OpenAI API usage that the repository should keep a separate Codex operations skill.

## Claude Code

- Anthropic's official overview positions Claude Code as a terminal-native coding tool and now highlights multiple Claude Code agents, the Agent SDK, GitHub Actions, GitLab CI/CD, channels, routines, and Chrome debugging.
- The official `anthropics/claude-code` repository currently deprecates npm installation in favor of install scripts, Homebrew, and Windows-native installers.
- The official `anthropics/claude-code-action` repository remains the main GitHub automation surface for Claude Code workflows.

## Gemini CLI

- Google's official `google-gemini/gemini-cli` repository presents Gemini CLI as an open-source terminal AI agent with MCP integration, Google Search grounding, checkpointing, headless scripting, GitHub integration, and `GEMINI.md` project context files.
- The releases section currently shows `v0.38.2` published on April 17, 2026, with explicit stable, preview, and nightly channels.
- This is a meaningful operational surface separate from ADK and generic Gemini API integrations.

## OpenAI Apps SDK

- The official Apps SDK docs define a complete ChatGPT app lifecycle centered on MCP apps, server setup, widgets, metadata, deployment, testing, submission, and reference documentation.
- The official Apps SDK home page states that approved apps appear in the ChatGPT apps store and that OpenAI creates a plugin for Codex distribution.
- The official examples repository is the strongest starting point for real widget plus MCP server composition.

## Implication for Repository Curation

- Provider APIs, agent frameworks, coding agents, and ChatGPT app surfaces should remain separate selection layers.
- Coding agents are now important enough to deserve first-class operational skills rather than scattered README mentions.
- Official-source reproducibility is improved when the repository keeps these surfaces indexed in both strategy docs and provider indexes.
