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
name: 'gram-mcp-gateway-operations'
description: 'Set up and troubleshoot Gram MCP gateways, register upstream servers, map enterprise identities to access rules, and verify policy enforcement with scoped access tests.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Gram MCP Gateway Operations

## Overview

Operate the MCP gateway in `speakeasy-api/gram`, now documented as Speakeasy Control Plane, from initial setup through access verification. Use the workflow to connect upstream tools with enterprise identity and scoped policies, then establish what the deployed configuration actually enforces rather than treating product security claims as guarantees.

## When to Use This Skill

- Set up a hosted organization or a local Gram development environment.
- Register a remote MCP server and expose it through a governed endpoint.
- Connect enterprise SSO and directory groups to server access rules.
- Configure and test MCP guardrails for a selected team or server.
- Investigate missing tools, authentication failures, unexpected access, or absent policy events.

## Core Capabilities

1. **Establish scope and install the appropriate environment.** Record the deployment URL, organization, project, release or commit, client, and intended change. Recheck the linked documentation against that deployment; feature availability can differ by plan and version. For hosted operation, use the organization setup board with `org:admin` access. For local development, use a checkout of the repository and its documented `./zero` bootstrap, which installs dependencies, migrates databases, and starts services; Mise and Docker are prerequisites. Inspect failures before retrying. This is a development bootstrap, not a verified production deployment recipe. Use `mise seed` only for disposable development data because rerunning it resets seeded data. See [repository setup](https://github.com/speakeasy-api/gram) and [organization onboarding](https://www.speakeasy.com/docs/ai-control-plane/getting-started).

2. **Register an upstream server.** In MCP Gateway > MCP, use Add new > Hosted remotely with a name and streamable HTTP URL; registration requires `mcp:write`. Verify connectivity checks reachability only: a 401 or 403 can pass. Inspect automatic OAuth setup results and configure missing credentials in the server's Identity settings. Keep caller authentication distinct from upstream credentials, and keep secrets out of reports. Use private availability for governed access; public servers skip server and per-tool access checks. Connect through Inspect to record tool annotations and review metadata changes before syncing. See [remote server operations](https://www.speakeasy.com/docs/ai-control-plane/mcp-gateway/remote-servers).

3. **Configure enterprise identity and role mapping.** In Organization settings > Team > IDP and SSO, verify the domain, configure SSO, then configure Directory Sync if entitled. Complete the provider portal steps and verify both connections. Map actual directory groups to intended roles and confirm a test user's resulting membership. Directory Sync makes provisioning and role assignment directory-controlled; do not attempt to repair it with manual invitations or manual member role assignment. Record a mapping of group identifier, platform role, intended server access, and test identity. See [IDP and SSO](https://www.speakeasy.com/docs/ai-control-plane/org-admin/identity).

4. **Apply server access rules.** Use Team Access with `org:admin`; examine Connect, View, and Manage together. Grants combine, while a block overrides grants for its capability. Unrestricted View or Manage can also satisfy connection checks. Remote servers and gateways support annotation-based narrowing, and unannotated tools can remain reachable; do not equate missing destructive annotations with safety. Review effective membership and conflicting roles before saving. Reload stale rules if a concurrent edit prevents the save. See [team access](https://www.speakeasy.com/docs/ai-control-plane/mcp-gateway/access/team-access).

5. **Assemble and connect a gateway.** Use Add new > New gateway, then Overview > Add servers. Copy the generated address into the client's supported MCP configuration and authenticate. Grant the intended identity access to both gateway and member servers: gateway access alone does not grant member access. Use live tool discovery rather than assuming a fixed gateway tool count. Check `connection_status` and finish each required upstream consent, then check again. A remote member's `unknown` health means unobserved, not necessarily failed. Signed-in gateways currently refuse agent API keys and workload sessions; verify the documented member URL or agent gateway route for unattended workloads. See [gateway endpoints](https://www.speakeasy.com/docs/ai-control-plane/mcp-gateway/gateway-endpoints).

6. **Configure MCP policy coverage explicitly.** Verify that MCP-scoped policies are enabled for the organization and that the operator has `org:admin`. Choose Selected MCP servers, detector, request/response surfaces, audience, and action. Select gateway endpoints explicitly; All MCP servers does not include them. Everywhere policies run through session instrumentation, not at the gateway. MCP policies support logging or denial, not warning confirmation or quarantine. Start with synthetic fixtures and logging, then test denial within the authorized scope. Document the published limit: scans inspect the first 50 KiB, and oversized content or scan timeouts are allowed through. Blocking a response cannot undo an upstream action already executed. See [guardrails and enforcement limits](https://www.speakeasy.com/docs/ai-control-plane/secure/guardrails).

7. **Verify access and diagnose failures.** Use separate allowed and disallowed test identities with fresh client sessions. Run the following acceptance checks against a harmless test server; never invoke a real destructive operation merely to prove denial.

   | Check | Evidence required |
   | --- | --- |
   | Allowed identity | Successful initialization, expected discovery, and one harmless tool result |
   | Disallowed identity | Restricted server or tool unavailable, plus rejection of a direct attempted invocation |
   | Policy match | Synthetic request or response produces the configured outcome and correlated risk event |
   | Policy control | Benign fixture succeeds without the policy falsely rejecting it |
   | Role removal | After directory convergence and reauthentication, former access is denied |

   Treat these as proposed acceptance tests, not claims that they have passed. Record endpoint, identity alias, tool, expected versus actual result, timestamp, and available event identifiers. Enable the required logging within the task's authorized retention and payload-capture settings. Inspect the member's tool logs when tracing routed calls.

   | Symptom | Next diagnostic step |
   | --- | --- |
   | Reachability passes but calls fail | Separate upstream credential failure from caller authentication; inspect OAuth metadata and consent |
   | Member or tool missing | Check gateway membership, availability, effective access, annotation restrictions, and client session freshness |
   | Unexpectedly broad access | Inspect public availability, overlapping grants, View/Manage permissions, and unannotated tools |
   | Policy does not block | Verify MCP scope, enabled state, audience, action, selected gateway, exclusions, size, and timeout behavior |
   | Identity role is stale | Check directory connection, group mapping, sync completion, and a fresh sign-in |
   | No evidence in logs | Confirm logging configuration, actual routed endpoint, time range, and member/gateway filters |

   If an access or policy test fails, stop distribution and report the failed invariant. Restore the recorded prior configuration when rollback is within scope; do not broaden permissions or disable policy merely to make a test pass.

## Inputs / Outputs

**Inputs:** Authorized deployment and project; deployment version; upstream URL and transport; credential references; IdP domain and group mappings; intended role/server/tool permissions; client configuration; policy intent; allowed and denied test accounts; harmless fixtures; prior settings for rollback.

**Outputs:** A redacted setup/change record containing server and gateway identifiers, credential mode, group-to-role mapping, access rules, policy scope/action, and the acceptance-test results. Include documentation URLs, observed limitations, unresolved prerequisites, and rollback status. Return client configuration only with credential placeholders.

**Completion criterion:** For an already provisioned environment, aim to complete one scoped configuration-and-verification cycle within 15 minutes. Success requires observed allowed access, denied unauthorized invocation, and the intended policy behavior with supporting evidence. If installation, entitlement activation, DNS verification, or directory provisioning remains pending, return a blocked or partial result with the exact next step; do not claim successful governance or a guaranteed rollout time.

## References

- [Source finding: speakeasy-api/gram](https://github.com/speakeasy-api/gram) — repository, product naming, and local development setup.
- [Getting started](https://www.speakeasy.com/docs/ai-control-plane/getting-started) — organization onboarding and logging prerequisites.
- [Remote MCP servers](https://www.speakeasy.com/docs/ai-control-plane/mcp-gateway/remote-servers) — registration, authentication, transport, and metadata.
- [IDP and SSO](https://www.speakeasy.com/docs/ai-control-plane/org-admin/identity) — enterprise sign-in and directory provisioning.
- [Team access](https://www.speakeasy.com/docs/ai-control-plane/mcp-gateway/access/team-access) — permission precedence and tool restrictions.
- [Gateway endpoints](https://www.speakeasy.com/docs/ai-control-plane/mcp-gateway/gateway-endpoints) — membership, client connections, and routed logging.
- [Guardrails](https://www.speakeasy.com/docs/ai-control-plane/secure/guardrails) — policy scope, supported actions, and enforcement limits.
