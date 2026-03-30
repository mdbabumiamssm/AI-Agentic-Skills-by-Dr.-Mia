<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA
-->

---
name: 'langsmith-observability-2026'
description: 'End-to-end tracing, evaluation, and incident response for agent workflows using LangSmith + OpenTelemetry.'
keywords:
  - observability
  - langsmith
  - tracing
  - open-telemetry
  - evals
measurable_outcome: 100% of mission-critical agents emit traces + evals with <5 min MTTR for regressions.
allowed-tools:
  - python
  - read_file
  - run_shell_command
---

# LangSmith Observability Skill

LangSmith unifies debugging, testing, deployment, monitoring, and evaluation for LLM+agent stacks. It now offers **end-to-end OpenTelemetry (OTel) support**, letting us standardize spans from LangChain/LangGraph runtimes and forward them to LangSmith or existing APM backends. [^otel]

IBM\'s 2026 overview summarizes why enterprises adopt LangSmith: consolidated debugging, detailed trace analysis, dataset-driven evals, and production monitoring. [^ibm]

## Integration Steps

1. **Instrument runtimes** – enable LangSmith tracing in `platform/core_kernel/workflow_engine.py` + `Agentic_AI/*` modules.
2. **OpenTelemetry bridge** – configure `LANGCHAIN_TRACING_V2=true` and `LANGSMITH_OTEL_EXPORTER=logging` (or OTLP) so traces reach both LangSmith + Splunk/Datadog.
3. **Eval pipelines** – define dataset YAML under `Skills/Research_Tools/LangSmith_Observability/evals/` mapping missions → metrics.
4. **Incident playbooks** – store runbooks in `ops/` describing how to triage high latency, hallucinations, or safety violations.
5. **Dashboarding** – pin LangSmith projects (Clinical Ops, Lab Automation, Research) and link to OpsGenie alerts.

## Sample Config (Python)

```python
import os
from langsmith import Client

os.environ.update({
    "LANGCHAIN_TRACING_V2": "true",
    "LANGCHAIN_PROJECT": "Clinical-Ops-Agent",
    "LANGCHAIN_API_KEY": "lsm-...",
    "LANGCHAIN_HUB_API_KEY": "lsm-hub-...",
    "LANGSMITH_OPENTELEMETRY_ENDPOINT": "https://otel.gateway",
})

client = Client()
trace_url = client.trace_url(run_id="{{run_id}}")
print(trace_url)
```

## Deliverables

| Path | Purpose |
|---|---|
| `dashboards/clinical_ops.json` | LangSmith snapshot exported for Grafana |
| `evals/acs_triage_evalset.yaml` | Baseline prompts + expected outputs |
| `ops/incident_runbook.md` | MTTR playbook for drift + safety regressions |

## References

[^otel]: LangChain blog (Mar 26, 2025) "Introducing End-to-End OpenTelemetry Support in LangSmith" describing bidirectional OTel integration for LangChain/LangGraph apps. 【turn0search13】
[^ibm]: IBM Think article (2026) "What is LangSmith?" outlining platform benefits (debugging, testing, monitoring) and enterprise readiness. 【turn2search11】

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
