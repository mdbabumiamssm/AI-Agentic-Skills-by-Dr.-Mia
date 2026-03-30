# LangSmith Observability Toolkit

## Folder Layout

```
LangSmith_Observability/
├── SKILL.md
├── dashboards/
├── evals/
├── ops/
└── scripts/
```

## TODOs

- [ ] `scripts/push_traces.py` – CLI for replaying stored JSON traces into LangSmith test projects.
- [ ] `dashboards/clinical_ops.json` – exported LangSmith dashboard with latency + hallucination metrics.
- [ ] `evals/acs_triage_evalset.yaml` – dataset feeding Regression Eval runs.
- [ ] `ops/incident_runbook.md` – response steps + contacts.

## Integration Tips

1. Keep LangSmith + OTel exporters enabled simultaneously.
2. Tag every run with `mission_id`, `patient_hash`, and `stack_version` for easier filtering.
3. Mirror critical traces to OpsGenie/PagerDuty for MTTA tracking.
