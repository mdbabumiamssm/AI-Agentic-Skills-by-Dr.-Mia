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
name: 'aws-agent-toolkit-integration'
description: 'Select, install, configure, and verify official AWS Agent Toolkit MCP servers, skills, or plugins for an AI coding agent working with AWS.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# AWS Agent Toolkit Integration

## Overview

Integrate an AI coding agent with the official AWS-supported Agent Toolkit by choosing the relevant plugin, local skill, or managed MCP connection. This workflow covers supported installation, credential and permission configuration, and one observable verification task without reproducing the upstream collection. Treat the 15-minute outcome as a workflow target with prerequisites available, not a measured performance guarantee.

## When to Use This Skill

- Connect an existing coding agent to AWS documentation or AWS APIs through the toolkit.
- Install task-specific AWS guidance or a plugin combining guidance with MCP configuration.
- Select a toolkit component for general AWS development, Bedrock/AgentCore, or analytics work.
- Diagnose an installation, authentication, permission, or component-discovery failure.
- Evaluate a migration from older AWS Labs integrations while preserving unrelated client settings.

## Core Capabilities

1. **Choose the smallest suitable component.** Read the current [toolkit README](https://github.com/aws/agent-toolkit-for-aws/blob/main/README.md) and the selected component documentation before changing configuration. Choose `aws-core` for general AWS development, `aws-agents` for Bedrock/AgentCore, or `aws-data-analytics` for S3 Tables, Glue, and Athena. Use standalone skills for local task guidance and a managed MCP connection for runtime tools. Skills and MCP operate independently; installing one does not prove the other is available.

2. **Use the supported installation path.** Identify the client, version, operating system, and intended project or user scope. Inspect existing configuration for duplicate integrations and preserve unrelated entries. Recheck upstream instructions because supported clients and commands can change. Examples documented by the toolkit:

   | Client or component | Installation route |
   | --- | --- |
   | Codex plugin | In a terminal, run `codex plugin marketplace add aws/agent-toolkit-for-aws`; then use `/plugins` inside Codex to install `aws-core`. |
   | Claude Code plugin | Inside Claude Code, run `/plugin install aws-core@claude-plugins-official`; substitute a documented plugin when needed. |
   | Cursor plugin | Import `aws/agent-toolkit-for-aws` through Settings → Plugins → Team Marketplaces → Add Marketplace → Import from Repo, then select the plugin. |
   | Standalone skills | Use `npx skills add aws/agent-toolkit-for-aws/skills`, selecting the intended client, skills, and installation scope. |
   | Direct MCP connection | Follow the linked AWS setup guide for the client's configuration format and supported authentication method. |

   Do not run slash commands in a shell. For a proxy connection, verify `uvx` is available; for skill installation, verify the runtime providing `npx`. Read each selected component's prerequisites rather than assuming the repository's Python language label defines every installation. Record the installed version or source revision when available. If Claude Code reports a missing plugin, refresh its official marketplace index with `/plugin marketplace update claude-plugins-official` and retry once.

3. **Configure authentication and permissions.** Follow the [AWS MCP setup guide](https://docs.aws.amazon.com/agent-toolkit/latest/userguide/getting-started-aws-mcp-server.html) for OAuth or a SigV4 proxy. OAuth suits supported remote clients; SigV4 supports local credential profiles and multi-account workflows. Configure the intended identity, endpoint, and operation Region; endpoint Region and default operation Region are distinct settings. For a local profile, verify identity with `aws sts get-caller-identity --profile <profile>` before account operations. Documentation search and local skill discovery do not require account credentials; authenticated API calls do.

   Use the organization's approved temporary credential or role workflow. Keep tokens, keys, and credential files out of prompts, checked-in configuration, and verification reports. Determine the permissions required by the selected operation and any connection-specific permissions from current documentation. Request only the needed actions and resources; never resolve an access denial by automatically granting administrator access. Apply documented IAM restrictions for agent actions when required by the task. Client tool visibility alone is not an IAM boundary.

4. **Verify a minimal task through the chosen component.** Reload or restart the client as its instructions require, then confirm the installed plugin, skill, or MCP connection appears. Discover actual tool names and schemas rather than guessing them. Choose one verification appropriate to the requested integration:

   | Integration | Minimal verification | Evidence required |
   | --- | --- | --- |
   | MCP documentation access | Ask the agent to find current AWS Lambda timeout documentation and retrieve the relevant page through its AWS documentation tool. | Successful tool invocation, official documentation URL, and a short answer supported by the retrieved content. |
   | Authenticated AWS API access | Through the connected MCP tool, execute an STS caller-identity lookup using the intended profile or identity. | Successful response containing `Account`, `Arn`, and `UserId`, with account and role matching the intended target. |
   | Standalone skill | Discover and load one installed skill, then apply its guidance to a small local planning task within its documented scope. | Skill identifier/path, evidence it loaded, and a concrete output following its instructions. |
   | Bundled plugin | Verify a bundled skill loads and perform one relevant MCP check above. | Both skill availability and successful runtime tool evidence. |

   A successful local AWS CLI identity check validates credentials, not the MCP integration. A documentation-only pass does not validate account access. Do not deploy resources or run security scans merely to prove installation; choose a task within the user's authorized scope. If the host cannot invoke MCP tools, report installation as configured but runtime verification as pending.

5. **Report results and recover narrowly.** Record the chosen component, installation scope, source/version, configuration location, authentication method without secrets, test input, observed output, and pass/fail status. Distinguish startup or transport errors from authentication errors, IAM denials, and incorrect Regions. Correct only the identified issue and retry the minimal task once; if it still fails, stop and report the blocker and the next concrete action. Preserve the prior configuration so a failed integration can be reversed without removing unrelated servers or skills.

## Inputs / Outputs

**Inputs:** Target agent/client and version; operating system; desired AWS task; installation scope; existing plugin/MCP/skill configuration; selected official component documentation; and, for authenticated work, the approved account, role/profile, Region, and permission constraints. Ask only for missing choices that affect installation or account access; never ask for secret values in chat.

**Outputs:** One selected integration with a reason for its selection; installed component or a precise prerequisite blocker; configuration changes with sensitive values excluded; and a verification report containing actual evidence. Mark each applicable check as passed, failed, or not run, and include elapsed time if measured. Completion requires successful execution through the selected component, not merely a configuration snippet or a proposed test.

## References

- [Agent Toolkit for AWS — official source repository](https://github.com/aws/agent-toolkit-for-aws).
- [Toolkit README — component selection and client installation](https://github.com/aws/agent-toolkit-for-aws/blob/main/README.md).
- [Official toolkit skills directory](https://github.com/aws/agent-toolkit-for-aws/tree/main/skills).
- [AWS MCP Server setup — authentication, configuration, and connection testing](https://docs.aws.amazon.com/agent-toolkit/latest/userguide/getting-started-aws-mcp-server.html).
- [AWS MCP Server tool reference](https://docs.aws.amazon.com/agent-toolkit/latest/userguide/understanding-mcp-server-tools.html).
