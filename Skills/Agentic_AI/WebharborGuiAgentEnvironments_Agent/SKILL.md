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



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'webharbor-gui-agent-environments'
description: 'Prepare WebHarbor website mirrors, reset task state, and record reproducible GUI agent evaluations using released Docker environments, Harbor tasks, and site verifiers.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# WebHarbor GUI Agent Environments

## Overview

This skill prepares controlled website environments, restores initial state, and evaluates GUI agents with traceable task and verifier artifacts. WebHarbor reports 153 Dockerized website environments and 3,518 tasks covering multi-step and authenticated workflows, addressing instability in live-web evaluation ([paper](http://arxiv.org/abs/2611.03634v1)). Use the released implementation to establish what is executable at the selected revision, and preserve evidence needed to reproduce each trial.

## When to Use This Skill

- Evaluate a GUI agent on a local WebHarbor mirror with repeatable initial conditions.
- Exercise authenticated workflows using the environment's seeded accounts and data.
- Diagnose differences between agent failures, verifier errors, and environment defects.
- Prepare isolated rollout environments for reinforcement learning.
- Extend a mirror for a new task while retaining existing workflow coverage.

## Core Capabilities

1. **Confirm and pin released artifacts.** Resolve environment, task, harness, and verifier versions before execution.
2. **Prepare executable mirrors.** Select a per-site Harbor container or the all-sites control-plane deployment.
3. **Restore task state.** Reset server state and establish a clean browser session before each independent trial.
4. **Specify auditable tasks.** Connect an instruction to its starting state, expected outcome, and released verifier.
5. **Record reproducible evaluations.** Preserve trajectories, state snapshots, verdicts, and configuration without exposing credentials.
6. **Evolve environments through tasks.** Validate functional changes against new and existing tasks, with human review of fidelity and task validity.

### Workflow

#### 1. Check availability and scope

The repository, root README, Harbor README, and startup script were accessible on 2026-11-03. They document public images and task datasets; this authoring check did not pull images or execute a benchmark. Recheck the actual artifacts needed for each run.

- Read the pinned checkout's `README.md`, `harbor/README.md`, and selected site's files under `sites/<site>/`.
- Record the repository commit, task dataset revision, task ID, image tag and resolved digest, and harness version. Mutable tags such as `latest` are insufficient to identify a reproducible run.
- Resolve the task through `ChilleD/WebHarbor-Tasks`: `registry.json` defines datasets and `tasks/` contains task directories. Inspect the selected directory before assembling commands.
- Confirm Docker availability, image access, free storage, browser readiness, and model endpoint reachability. The documented all-sites image is about 20 GB; cold downloads or builds can exceed the time budget.
- If an image, task, verifier, or required asset is unavailable, report the exact missing artifact and stop dependent execution. Do not invent substitute interfaces or claim a completed evaluation.

The 15-minute target covers one prepared task, reset validation, and an evidence report with cached dependencies and an available model endpoint. Full-suite evaluation, cold setup, and RL training are separate jobs. If the budget expires, return an explicitly incomplete report with the last verified stage.

#### 2. Select the deployment mode

**Harbor task mode:** Prefer this for graded evaluations. Released tasks use `wildawm/webharbor:<site>-v1`; per-site images serve their registered port and ship without verifiers. The task package supplies verification code. Install `harbor==0.24.0` for the documented harnesses; the Harbor README identifies an API incompatibility with 0.25.

Run from the checkout's `harbor/` directory with `PYTHONPATH=agents`. Inspect `harbor run --help` for the installed version and select exactly the intended task using its supported local-task or filtering option. Do not launch the whole test split for a single-task smoke test.

For the paper's GUI harness, use import path `webharbor_agents.openwebrl:OpenWebRLWebHarbor` with `--ak api_mode=chat_completions`, an appropriate `-m openai/<model>` value, and the documented `OPENAI_BASE_URL` and `OPENAI_API_KEY` environment variables. Resolve placeholders from the user's configuration. The endpoint must be reachable inside the container; container `localhost` does not address a model server on the Docker host.

**All-sites mode:** Use for interactive inspection or sequential rollouts that need the documented reset service. The following example adapts the released quickstart to bind host ports to loopback:

```bash
export WEBSYN_CONTROL_TOKEN="$(python3 -c 'import secrets; print(secrets.token_urlsafe(48))')"
docker run --name webharbor-eval -d \
  -e WEBSYN_CONTROL_TOKEN \
  -p 127.0.0.1:8101:8101 \
  -p 127.0.0.1:40000-40152:40000-40152 \
  wildawm/webharbor:all-v1
```

Use a unique container name and record its image digest. `WEBSYN_CONTROL_TOKEN` must contain at least 32 characters; keep it out of agent prompts and saved reports. Port assignment follows the `SITES` array in `websyn_start.sh`, starting at 40000; for example, Amazon is on 40001. Verify the selected site's response and startup logs before beginning; the control plane starts on 8101 only after all sites are ready.

For source builds, the documented entry point is `./scripts/build.sh webharbor:dev`, which retrieves missing static assets from `ChilleD/WebHarbor`. Record asset revisions and build results as well as the source commit.

#### 3. Establish and verify initial state

For an all-sites Amazon trial, the released reset interface is:

```bash
curl --fail-with-body -sS -X POST \
  -H "Authorization: Bearer $WEBSYN_CONTROL_TOKEN" \
  http://localhost:8101/reset/amazon
```

Use the registered site identifier for other sites. `/reset-all` exists, but use it only when the entire deployment belongs to this run and no other rollout is active. Do not assume per-site Harbor containers expose this same control plane; inspect their task lifecycle and entrypoint, and let Harbor create each task environment.

After reset, verify an observable seed-state invariant from the selected task, such as the expected cart contents or account records. Clear prior browser cookies, storage, and session state, then establish the task's documented authentication precondition using seeded local credentials. Server reset alone does not prove browser isolation.

Before trusting repeated rollouts, perform one reversible local mutation, reset, and check that the seed invariant returns. Save this result. Never reset a shared mutable site during another trial; allocate independent containers for concurrent work.

#### 4. Bind the task to its evaluator

Read the released task instruction, environment configuration, `tests/test.sh`, and relevant site task entry before execution. Record the task ID, local start URL, authentication precondition, expected persistent changes or final answer, action/time budget, and verifier identity. Keep grading code and expected answers outside the agent's observation surface.

For new tasks, follow the actual schema of a neighboring released task rather than inventing a WebHarbor schema. Add concrete success conditions and a verifier covering meaningful state changes or evidence. Validate that the task is achievable from the reset state, and check that an incomplete action does not receive a successful verdict.

Keep training and held-out task selection explicit. Use the registry's published split membership; do not silently move held-out tasks into training or assume every split contains unseen websites.

#### 5. Execute and preserve evidence

Run the selected GUI harness with fixed task selection, observation mode, model revision, sampling configuration, and budgets. Record any supported random seed; do not claim that deterministic website state makes a stochastic model deterministic.

In the documented Harbor verification flow, `tests/build_trajectory.py` reconstructs the run evidence, including `initial.db`, `after.db`, access-log steps, screenshots, final answer, and termination status. `tests/test.sh` invokes the site's `verify_<n>.py` and writes `/logs/verifier/reward.json`, containing binary reward and pass/fail check counts. Preserve the resulting trial artifacts before container disposal.

Inspect the reward and verifier logs together. A missing or malformed verdict is an evaluation error, not proof of agent failure. Report environment or task defects separately, retain failed attempts, and disclose any retries. If using the optional agent judge, retain its verdict separately: the released judge audits but does not replace the official verifier reward.

Repeat a selected task from a fresh reset when checking reproducibility. Compare initial-state invariants and evaluator behavior; report model-output variation without attributing it automatically to environment drift.

#### 6. Extend and review when requested

For a new workflow, identify the missing UI, backend operation, or seed data, implement the smallest supporting change, and rerun both the new task and affected existing tasks. Request human review of visual fidelity, functional correctness, and task validity before treating the extension as validated. Record a new environment/task revision and keep evaluation results from different revisions distinguishable.

## Inputs / Outputs

**Inputs:** A repository revision; site and task identifiers; pinned task dataset and image; Docker runtime; model endpoint and secret references; harness choice; authentication preconditions; evaluation budget; and an output location. A prepared single-task run is the default scope.

**Outputs:** A run manifest, reset-check evidence, task specification, trajectory and screenshots, initial/final state evidence, raw verifier verdict and logs, and a concise evaluation report. These are workflow deliverables, not claims that WebHarbor defines a particular manifest schema.

The manifest should capture:

- Run ID and timestamps; code commit, image digest, dataset revision, task ID, and split.
- Browser/harness versions, model identity, observation mode, sampling settings, budgets, and supported seed.
- Local site address, reset method, seed-state checks, and authentication setup without credentials.
- Exact invocation with secrets redacted, artifact locations, reward/check counts, termination status, and error classification.

A valid completed output identifies the selected task and artifact versions, demonstrates readiness and initial-state restoration, includes an interpretable verifier result with supporting evidence, and states whether the trial passed. Setup failure, timeout, or absent evidence must remain explicitly incomplete. For multi-task runs, report attempted, completed, failed, and errored counts and the denominator used for any success rate.

## References

- [WebHarbor: Docking Real Websites for Evolving GUI Agent Environments, arXiv:2611.03634v1](http://arxiv.org/abs/2611.03634v1) — source finding and reported scope.
- [WebHarbor repository and quickstart](https://github.com/aiming-lab/WebHarbor/blob/main/README.md) — images, build instructions, reset routes, and port lookup.
- [WebHarbor Harbor integration](https://github.com/aiming-lab/WebHarbor/blob/main/harbor/README.md) — task packaging, harness compatibility, and verifier evidence.
- [WebHarbor startup implementation](https://github.com/aiming-lab/WebHarbor/blob/main/websyn_start.sh) — seed restoration, registered site order, and readiness checks.
- [Released Harbor tasks](https://huggingface.co/datasets/ChilleD/WebHarbor-Tasks) — dataset linked by the repository.
- [Environment assets](https://huggingface.co/datasets/ChilleD/WebHarbor) and [Docker images](https://hub.docker.com/r/wildawm/webharbor) — distribution endpoints linked by the repository; verify access at execution time.
