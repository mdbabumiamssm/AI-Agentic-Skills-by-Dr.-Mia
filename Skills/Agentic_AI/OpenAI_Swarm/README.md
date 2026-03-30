# OpenAI Swarm Playbook (2026)

| Component | Purpose | File |
|---|---|---|
| Skill definition | Usage, outcomes, guardrails | `SKILL.md` |
| Blueprint templates | Reusable agent instructions & tools | `playbooks/swarm_blueprints.md` |
| Demo runner | CLI entry point that wires Swarm into in-repo toolchains | `scripts/run_swarm_demo.py` |

## Architecture Snapshot

1. **Agents** – YAML/JSON definitions with `instructions`, `functions`, `context_variables`.
2. **Handoffs** – Tools return `Agent` objects. Each handoff is auditable for clinical workflows.
3. **Routines** – Pre-packaged sequences (triage → labs → safety) triggered by `swarm.routines`.
4. **Upgrade Path** – Export stable Swarm missions to OpenAI Agents SDK using shared instruction files.

## Backlog

- [ ] Add `swarm_blueprints.md` with templated agent JSON for hematology, oncology, and self-driving lab orchestration.
- [ ] Implement `scripts/run_swarm_demo.py` bridging Swarm with `platform/adapters/runtime_adapter.py`.
- [ ] Wire Swarm traces into LangSmith via OTel exporters (see `Research_Tools/LangSmith_Observability`).
