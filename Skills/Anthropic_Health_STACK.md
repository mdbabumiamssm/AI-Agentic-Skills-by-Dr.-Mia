# Anthropic Health & Life Sciences Stack (2026)

**Focus:** Operationalizing Anthropic's \"Briefing: Healthcare and Life Sciences\" vision inside the Universal Biomedical Skills Platform. This stack emphasises asynchronous coworker agents, MCP-native tooling, and regulatory-grade audit trails.

## Components

| Layer | File(s) | What It Does |
| --- | --- | --- |
| Inbox Router | `Clinical/anthropic_inbox_router.py` | Parses work requests (prior auth, safety, regulatory) and fans them out to specialist coworkers via the Event Bus. |
| Prior Authorization Coworker | `Clinical/Prior_Authorization/anthropic_coworker.py` | Emits `<thinking>`, `<analysis>`, `<decision>` blocks plus structured JSON justifications to mirror Anthropic coworker transcripts. |
| Regulatory Response Coworker | `Pharma/Regulatory_Affairs/anthropic_regulatory_coworker.py` | Drafts and reviews CTD/Module 2 updates with built-in policy citations and routing for human sign-off. |
| Pharmacovigilance Monitor | `Clinical/Safety/pharmacovigilance_monitor.py` | Prioritizes adverse event signals, tags risk levels, and publishes asynchronous status cards. |

## Operating Pattern

1. **Event Intake** – The inbox router listens to `anthropic.health.inbox` topics (e.g., \"PA request\", \"safety signal\", \"FDA query\").
2. **Coworker Loop** – Specialist modules fetch artifacts, run reasoning traces, emit `<thinking>` logs, and deposit deliverables back on the bus.
3. **Audit & Replay** – Every coworker attaches citations, policy versions, and timestamped determinations, enabling full audit replay.
4. **Human in the Loop** – Each worker can request clarification or escalate via `needs-human-review` events, matching Anthropic's \"copilot/coordinator\" messaging.

## Usage

```bash
# Kick off a prior authorization coworker run
python3 Skills/Clinical/Prior_Authorization/anthropic_coworker.py examples/mri_case.json

# Route multiple inbox items (mock queue)
python3 Skills/Clinical/anthropic_inbox_router.py examples/inbox.json
```

Each module prints the structured payload that would be sent to Claude (XML style) along with the JSON representation that flows through BioKernel.

## Integration Notes

* Register MCP tools (`anthropic.prior_auth`, `anthropic.regulatory`, `anthropic.pv_monitor`) so Anthropic desktop clients can subscribe directly.
* Extend the Optimizer with an `anthropic` target that wraps tasks inside the `<thinking>/<analysis>/<decision>` contract.
* Mirror Anthropic's compliance stance by persisting every coworker exchange to the Audit Log and signing responses with the Event Bus digest helper.
* Drive every coworker definition from USDL specs so the Anthropic prompt/tag layout always matches the OpenAI/Gemini variants (see `docs/USDL_OVERVIEW.md`).

This stack keeps us aligned with Anthropic's emphasis on safe coworker agents that respect regulatory guardrails while still moving fast on administrative heavy-lifting. 
