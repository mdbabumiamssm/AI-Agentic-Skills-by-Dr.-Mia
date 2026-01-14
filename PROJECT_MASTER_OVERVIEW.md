# Universal Biomedical Skills Platform - Master Project Overview

**Last Updated:** January 14, 2026
**Status:** Active Development (Enterprise Edition v2026.3)

## 1. Project Identity & Mission
This project, the **Universal Biomedical Skills Platform**, is a comprehensive "Agentic Operating System" for Life Sciences. It standardizes how AI agents (running on OpenAI, Anthropic, Gemini, or open models) collaborate to solve complex scientific problems.

**Core Philosophy:** "Write Once, Run Anywhere" for biomedical AI skills, capable of orchestration via a "BioKernel".

## 2. Key Components & Architecture

### A. The Platform (`/platform`)
The core runtime environment ("BioKernel") that executes skills.
-   **BioKernel (`platform/biokernel`):** A FastAPI/MCP-based engine that manages agent orchestration and event buses.
-   **Optimizer (`platform/optimizer`):** Features a "USDL Transpiler" (Universal Skill Definition Language) that converts generic skill specifications into provider-optimized prompts (e.g., specific JSON schemas for OpenAI or XML thinking blocks for Anthropic).
-   **Dashboard (`platform/dashboard.py`):** Real-time visualization of agent workflows.
-   **Adapters (`platform/adapters`):** Connectors for Claude, OpenAI, and potentially Gemini.

### B. The Skills Library (`/Skills`)
A massive collection of domain-specific capabilities, organized by scientific discipline.
-   **Dual-Stack Alignment:**
    -   **OpenAI Health Stack:** JSON Schema-enforced agents for clinical ops and triage.
    -   **Anthropic Health Stack:** "Thinking Block" agents for high-stakes regulatory/safety tasks.
-   **Domains:** Genomics, Drug Discovery, Clinical Trials, Immunology, etc.
-   **Agentic AI:** Contains the `SwarmOrchestrator` for multi-agent teams.

### C. Skill Collections (`/skill collections`)
External repositories and reference implementations (e.g., Auto-GPT, LangChain, MCP Servers) used for inspiration or direct integration.

## 3. Current Strategic Focus (Jan 2026)
1.  **Antigravity Integration:** Adopting Google's "Antigravity" platform standards (released Nov 2025), specifically the `SKILL.md` format and native MCP integration.
2.  **MCP (Model Context Protocol):** Deepening support for MCP to allow agents to interact directly with external tools and data contexts standardly.
3.  **Swarm Orchestration:** Moving from single-agent tasks to multi-agent "swarms" (Researcher, Reviewer, Safety Officer).

## 4. How to Navigate
-   **To run the kernel:** `python platform/biokernel/server.py` (check specific instructions).
-   **To see the roadmap:** Read `SUBSTANTIAL_IMPROVEMENT_PLAN_2026_UPDATED.md`.
-   **To add a new skill:**
    1.  Define it in `Skills/<Domain>/`.
    2.  Use the `USDL` format or the new `SKILL.md` format.
    3.  Register it with the BioKernel.

## 5. Maintenance
*This document should be updated whenever major architectural changes occur or new top-level directories are added.*

## 6. Antigravity & MCP Integration (New)
**Released Jan 2026** - We have integrated Google Antigravity standards.

### A. Antigravity Skills (`skill collections/Antigravity_Skills`)
We now support the `SKILL.md` format.
-   **Structure:** YAML frontmatter (`name`, `description`) + Markdown body (`Instructions`).
-   **Current Skills:**
    -   `clinical-trial-matcher`: Matches patients to trials.
    -   `crispr-designer`: Designs gRNAs with off-target analysis.
    -   `regulatory-drafter`: Anthropic-style regulatory writing.

### B. MCP Server (`platform/biokernel/mcp_server.py`)
A JSON-RPC 2.0 compliant server compatible with Claude Desktop and other MCP clients.
-   **Tools:** Exposes `run_bio_agent` to execute any registered skill.
-   **Prompts:** Exposes all Antigravity skills as system prompts.
-   **Transport:** Standard Input/Output (stdio).

**To run the MCP Server:**
```bash
python3 platform/biokernel/mcp_server.py
```
*(Ensure `fastapi` and `uvicorn` are installed)*
