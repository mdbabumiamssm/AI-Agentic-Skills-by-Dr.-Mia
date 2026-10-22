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
name: 'mcp-gateway-registry-operations'
description: 'Deploy and operate MCP Gateway & Registry, register servers, configure Keycloak or Entra OAuth identities, and verify scoped tool access and audit records.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# MCP Gateway and Registry Operations

## Overview

Operate the MCP Gateway & Registry to centralize MCP server discovery and authenticated access for agents and coding assistants. This skill connects deployment, identity configuration, server registration, and fine-grained authorization with observable access and audit checks. Use documentation matching the deployed revision because API schemas, scope storage, and authentication behavior can change between releases.

## When to Use This Skill

- Deploy or validate an instance using Docker Compose, AWS ECS Terraform, or Amazon EKS Helm charts.
- Register an MCP server and make its tools discoverable through the gateway.
- Configure Keycloak or Microsoft Entra ID for human users or machine identities.
- Restrict a principal to selected servers, MCP methods, and named tools.
- Diagnose authentication failures, unexpected permissions, missing tools, or absent audit events.
- Verify an existing configuration after a deployment, identity, or scope change.

## Core Capabilities

1. **Select and validate a supported deployment.**
   Record the repository tag or commit, deployment target, gateway URL, storage backend, and current configuration before changing anything. Read the matching installation and configuration guides; use the ECS Terraform or EKS Helm documentation only for that target. The documented Compose quick start uses `.env.example` and `./build_and_run.sh --prebuilt`; configure required secrets and initialize MongoDB/DocumentDB and identity services using that revision's guide first. Preserve an existing environment file and use the documented health checks to verify the registry, auth server, proxy, and database independently.

2. **Configure the identity provider and authentication flow.**
   Use `AUTH_PROVIDER=keycloak` or `AUTH_PROVIDER=entra` as documented. For Keycloak, verify realm, clients, redirect URIs, group claims, and the service-account configuration when M2M is required. For Entra, follow its setup guide for tenant/application configuration, redirect URIs, audience, consent, and group identifiers; do not substitute group display names where object IDs are required. Distinguish browser authorization-code sessions from workload client-credentials tokens and the project's programmatic-token flow. Keep client secrets, bearer tokens, and cookies out of reports and shell tracing; validate issuer, audience, expiry, and signature through the configured authentication path.

3. **Register and discover MCP servers through documented interfaces.**
   Inspect `api/openapi.json`, the API guide, and the selected revision's client or CLI before constructing a request. First list existing registrations to avoid creating a duplicate path. Use the Registry UI, `api/registry_client.py`, `api/registry_management.py`, or the documented registration API; preserve unrelated records.

   The API guide documents these interfaces; confirm their request encoding and authentication against the installed version:

   | Operation | Documented interface | Verification |
   |---|---|---|
   | List registrations | `GET /api/servers` | Locate the intended server and check for path collisions. |
   | Register | `POST /api/servers/register` | Form fields include `name`, `description`, `path`, and `proxy_pass_url`; documented success is HTTP 201. |
   | Read registration | `GET /api/servers/{path}` | Confirm the stored route and backend URL. |
   | Inspect discovered tools | `GET /api/tools/{service_path}` | Compare discovered names with the backend's expected tools. |
   | Refresh health and tools | `POST /api/servers/{path}/rescan` | Inspect health and discovery results after the refresh. |
   | Enable or disable | `POST /api/servers/toggle` | Form fields are `path` and `new_state`; read back the resulting state. |

   Choose a backend URL reachable from the gateway's network; container-local `localhost` does not identify a separate container. Verify the registered transport and endpoint path. Configure upstream credentials separately from gateway login, using the documented credential mechanism. Registration success alone does not prove that scanning, activation, discovery, or invocation succeeded; inspect each relevant state before testing access.

4. **Apply scoped access without conflating management and execution.**
   Resolve group-to-scope mappings for human users and the documented scope assignments for workload identities. Registry UI/API permissions and MCP server permissions are separate. Grant only the requested server, the protocol methods needed by the client, and the named tools needed by the task. Calling a tool requires both `tools/call` permission and authorization for that tool; listing tools does not prove execution permission.

   Follow the scope-management API or documented seed/import workflow. Current documentation stores scope documents in the MongoDB/DocumentDB `mcp_scopes` collection and notes that legacy `scopes.yml` was removed in v1.24.8; do not treat YAML illustrations as a live configuration file. Check all effective scopes for broader grants before concluding that a restricted scope enforces the intended boundary. If virtual servers are in scope, check the selected revision's virtual-server rules, including backing-server grants and alias/original tool-name matching where required.

5. **Verify access using the actual client identity.**
   Use an MCP-capable client with the configured transport to perform initialization, send `notifications/initialized`, list tools, and call one authorized, harmless tool. Preserve session headers required by that transport. Test through the gateway URL, and separately diagnose backend reachability if needed.

   Run this small verification matrix using designated test identities and safe arguments:

   | Case | Expected observation |
   |---|---|
   | Intended principal, allowed tool | Successful MCP result with the expected benign output. |
   | Same principal, excluded tool | Authorization denial; no successful backend action. |
   | Principal without the server grant | Server access denied. |
   | Missing or invalid credential | Authentication rejected before successful tool execution. |

   Check HTTP and MCP JSON-RPC outcomes: HTTP 200 alone is not proof of tool success. Record the actual status and error instead of assuming every denial uses the same code. Do not use an administrator token to validate restricted-user behavior. If an excluded operation succeeds, stop further calls, inspect effective grants and routing, and correct only the authorized configuration.

6. **Verify audit evidence and produce an operational result.**
   Use the admin audit UI or documented read interfaces, including `GET /api/audit/events?stream=registry_api&limit=50` and `GET /api/audit/events?stream=mcp_access&limit=50`. Match the verification window, identity, server, method/tool, request identifiers where available, and reported outcomes. Check gateway/auth logs for rejected requests when the selected audit stream does not capture that rejection; state the coverage gap explicitly.

   Confirm that credentials are masked and that the configured durable sink receives new records. Inspect audit enablement, runtime write errors, retention settings, and initialized TTL indexes if events are missing or unexpectedly retained. Do not infer complete delivery, tamper evidence, or compliance from audit enablement alone. Preserve existing retention unless changing it is part of the task; shortening retention can delete older records.

   For failures, distinguish identity validation, scope mapping, proxy routing, backend authentication, discovery/scanning, and audit storage before applying a fix. Read back changed records and restore the prior configuration if a failed change must be rolled back within the authorized scope. Report a missing prerequisite or unresolved failure rather than weakening authorization to make a check pass.

## Inputs / Outputs

**Inputs**

- Requested operation, authorized environment, repository revision, and deployment mode.
- Gateway URL, network/TLS configuration, and MongoDB or DocumentDB configuration.
- Identity-provider details, secret references, intended groups/scopes, and test principals.
- MCP server name, registration path, backend URL, transport, upstream-auth requirements, and expected tools.
- One harmless allowed invocation, one excluded operation to probe, and access to the relevant audit view.

**Outputs**

- A deployment or configuration result tied to the inspected revision, with changes and unresolved prerequisites identified.
- Registration and discovery evidence, plus the effective principal-to-server/method/tool access mapping.
- A pass/fail/not-run verification matrix with redacted request details, HTTP/MCP outcomes, and audit references.
- Any required rollback details and a concise next action for unresolved failures; never include credentials.

**Completion criterion:** On a provisioned environment, complete the requested operation and its focused verification within the 15-minute target, with read-back evidence, an allowed call, applicable denial checks, and audit evidence. Initial infrastructure provisioning or identity approvals may exceed that target; report incomplete checks as not run and do not claim successful deployment from a plan alone.

## References

- [Source repository: MCP Gateway & Registry](https://github.com/agentic-community/mcp-gateway-registry).
- [Installation](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/installation.md) and [configuration](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/configuration.md).
- [AWS ECS Terraform deployment](https://github.com/agentic-community/mcp-gateway-registry/blob/main/terraform/aws-ecs/README.md) and [EKS Helm charts](https://github.com/agentic-community/mcp-gateway-registry/blob/main/charts/README.md).
- [Authentication](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/auth.md), [complete setup](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/complete-setup-guide.md), and [Entra ID setup](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/entra-id-setup.md).
- [API guide](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/registry_api.md) and [OpenAPI specification](https://github.com/agentic-community/mcp-gateway-registry/blob/main/api/openapi.json).
- [Access control and scopes](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/scopes.md) and [scope management](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/scopes-mgmt.md).
- [Dynamic tool discovery](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/dynamic-tool-discovery.md) and [coding-assistant integration](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/ai-coding-assistants-setup.md).
- [Audit logging and its limitations](https://github.com/agentic-community/mcp-gateway-registry/blob/main/docs/audit-logging.md).

Reference links target the upstream default branch. Resolve them at the deployed tag or commit before executing version-sensitive operations.
