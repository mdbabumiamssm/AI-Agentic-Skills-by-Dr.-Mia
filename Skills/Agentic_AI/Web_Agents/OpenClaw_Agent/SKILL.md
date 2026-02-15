---
name: OpenClaw_AI_Agent
description: Autonomous local AI agent ("The Personal AI OS") capable of controlling browsers, sending emails, and executing complex workflows via TypeScript/Python bridges.
version: 1.0.0
category: Agentic_AI
author: AI Agentic Skills Team
input:
  task: Natural language command (e.g., "Find the cheapest flight to Tokyo and email me the summary").
  mode: "headless" or "ui" (interactive).
output:
  result: Execution log or summary.
  artifacts: Generated files or screenshots.
---

# OpenClaw AI Agent

**OpenClaw AI** (formerly Clawdbot/Moltbot) is a conversation-first, autonomous agent that runs locally on your machine. It bridges the gap between LLMs and your local operating system, allowing for secure, private, and powerful automation.

## Capabilities
1.  **Browser Control:** Uses Playwright/Puppeteer protocols to navigate websites, click buttons, and scrape dynamic content.
2.  **System Integration:** Can write files, execute shell commands, and interact with local applications.
3.  **Vision Support:** Analyzes UI elements using vision-capable models (e.g., GPT-4o, Claude 3.5 Sonnet) to "see" the screen.
4.  **Privacy First:** All data stays local (SQLite) unless explicitly sent to an external LLM provider.

## Usage (Python Wrapper)

```python
from Skills.Agentic_AI.Web_Agents.OpenClaw_Agent.openclaw_wrapper import OpenClaw

claw = OpenClaw(headless=False)
result = claw.run_task("Go to github.com/trending and list the top 3 python repos.")
print(result)
```

## Security Note
⚠️ **Warning:** OpenClaw executes real commands on your machine. Ensure you run it in a sandboxed environment (Docker/VM) if testing untrusted prompts.
