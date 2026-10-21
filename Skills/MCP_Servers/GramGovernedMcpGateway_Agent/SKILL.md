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
name: 'gram-governed-mcp-gateway'
description: 'Configure and validate Gram MCP gateway access through upstream registration, identity-provider role mapping, and per-action policy checks, with deployment-specific feature verification.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Gram Governed MCP Gateway

## Overview

Use this skill to establish and verify a governed path from agents through Gram to upstream MCP servers. Gram describes role-based distribution through enterprise identity providers and policy enforcement on agent actions; this workflow turns those capabilities into a scoped configuration and observable access checks. Verify the selected deployment and its entitlements before treating advertised capabilities as available.

## When to Use This Skill

- Set up or assess a Gram gateway for a team or project.
- Register an upstream MCP server and expose selected tools to agents.
- Map identity-provider groups or roles to gateway access.
- Configure tool access or action policies and compare allowed and denied calls.
- Investigate unexpected tool visibility, authorization failures, or policy bypasses.
- Determine whether a requested control is available in an open-source deployment or requires a separately licensed service or component.

## Core Capabilities

1. **Deployment and entitlement discovery.** Record the version or commit, deployment mode, project boundary, and available administration interface. The repository root uses AGPL-3.0; its README describes an open-source stack and a hosted platform. These facts do not establish a complete feature-to-license matrix. Check component licenses, version-specific documentation, and deployment entitlements; mark each requested control as verified available, separately licensed with evidence, unavailable, or unverified. Do not assume IdP synchronization or policy features are exclusively paid or universally included.

2. **Gateway setup.** Follow documentation matching the selected release. For local development, the README identifies `./zero` as the bootstrap entry point, with Mise and Docker dependencies; inspect the script before executing it in an authorized checkout. This is a development workflow, not evidence of a production deployment recipe. Reuse an existing instance when supplied and confirm readiness before changing access configuration.

3. **Upstream registration.** Capture the MCP endpoint, supported transport, authentication method, and intended tool subset. Register it through the documented interface for the running version, reference credentials through the deployment's secret mechanism, and check connectivity and tool discovery. Keep client-to-gateway authentication distinct from gateway-to-upstream credentials. Never invent API routes, CLI flags, or configuration keys.

4. **Identity mapping.** Identify the IdP, trusted issuer, relevant claims, and intended group-to-role mapping using the supported integration. Prepare at least one permitted identity and one authenticated identity without the required role. Verify effective membership and project scope from the gateway rather than assuming that a successful login grants tool access. If integration is unavailable, report the gap; do not substitute shared administrator credentials as proof of role provisioning.

5. **Access and action policy configuration.** Translate the requested access matrix into supported team, server, and tool permissions. Configure only documented policy mechanisms and scope each rule to the intended identities and actions. Keep a recoverable record of the previous configuration. An advertised inspection capability is not proof that a specific blocking rule is enabled or enforced.

6. **Allowed-versus-denied validation.** Run the following checks through the gateway using harmless test tools or approved fixtures. Record identity, tool, expected result, observed result, timestamp, and available request or audit identifiers.

   | Check | Expected evidence |
   | --- | --- |
   | Permitted identity calls an allowed tool | Successful upstream result and attributable gateway evidence |
   | Authenticated identity without access calls the same tool directly by name | Authorization denial before upstream execution |
   | Permitted identity calls a forbidden tool | Denial even if the caller knows the tool name |
   | Missing or invalid authentication | Rejection under the configured authentication requirement |
   | Synthetic input matches an enabled blocking policy | Policy-specific denial with no upstream side effect |

   Compare `tools/list` visibility with `tools/call` enforcement; hiding a tool alone is insufficient. Distinguish authorization denial from connection failure, upstream errors, or invalid arguments. Use upstream logs or fixture counters where available to establish that rejected calls did not execute; otherwise label that evidence incomplete. Mark unsupported policy tests as not run, not passed.

7. **Evidence and recovery.** Summarize configuration changes, verified controls, failed checks, and entitlement gaps. If an expected denial succeeds, stop further exposure and restore the affected configuration within the authorized scope. Do not weaken authentication or broaden roles to make a test pass. Redact tokens and sensitive payloads from all output.

## Inputs / Outputs

**Inputs:** Deployment URL or authorized checkout; version; administrative access; upstream endpoint and credential references; IdP configuration and test identities; desired role/server/tool access matrix; policy requirements; harmless test fixtures; and the permitted change scope.

**Outputs:** A redacted configuration summary, upstream registration status, effective identity mapping, feature/entitlement evidence, and a validation matrix with pass, fail, or not-run results. Include unresolved prerequisites and restoration instructions for changes made.

**Execution boundary:** Target a 15-minute configuration and verification pass when the instance, credentials, and test identities are ready. Fresh infrastructure or IdP provisioning may exceed that window; report a concrete prerequisite gap instead of claiming completion. A successful governance result requires the allowed call and applicable denied-call checks to behave as expected, with sufficient evidence to distinguish enforcement from unrelated failures.

## References

- [Gram repository and source finding](https://github.com/speakeasy-api/gram) — product scope and entry point for current documentation.
- [Gram README](https://github.com/speakeasy-api/gram/blob/main/README.md) — platform capabilities and local development bootstrap.
- [Gram root license](https://github.com/speakeasy-api/gram/blob/main/LICENSE) — repository licensing; inspect component-specific notices and deployment terms separately.
- [Gram documentation source](https://github.com/speakeasy-api/gram/tree/main/docs) — locate version-relevant configuration instructions before execution.
