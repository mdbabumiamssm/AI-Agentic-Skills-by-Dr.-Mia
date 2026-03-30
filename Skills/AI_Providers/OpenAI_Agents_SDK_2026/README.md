# OpenAI Agents SDK Integration Plan

| Track | Owner | Deliverable |
|---|---|---|
| SDK parity | Platform Core | `notebooks/openai_agents_sdk_quickstarts.ipynb` |
| Guardrails | Safety & Compliance | `platform/adapters/security_monitor.py` hooks for AgentKit Guardrails |
| UI | Product | ChatKit embed + mission control page |
| Migration | Docs | `docs/Assistants_v1_to_AgentsSDK.md` |

## Milestones

1. **M0 – Sandbox**: run local Agents SDK sample, ensure `Responses API` + `web_search` + `file_search` works with MCP registry.
2. **M1 – Mission Templates**: export templates for Clinical Ops, Lab Automation, and Research Agents.
3. **M2 – Observability**: connect traces/evals to `Research_Tools/LangSmith_Observability` dashboards.
4. **M3 – Enterprise Rollout**: enable Connector Registry with approved data endpoints + Guardrails.
