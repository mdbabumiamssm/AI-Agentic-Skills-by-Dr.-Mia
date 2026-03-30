# LangSmith Incident Runbook (Draft)

1. **Detect** – PagerDuty alert triggered by LangSmith eval regression or latency > p95.
2. **Triage** – Open LangSmith project, filter by `mission_id`, inspect traces + spans.
3. **Contain** – Toggle `SAFE_MODE=true` to disable autonomous execution, fall back to manual review.
4. **Diagnose** – Compare eval history vs. new run; check tool latency vs. hallucination tags.
5. **Fix** – Roll back stack_version or deploy patch; rerun evalset.
6. **Review** – File RCA in `docs/ops/postmortems/<date>.md`.
