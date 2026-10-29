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
name: 'hol-guard-agent-runtime-security'
description: 'Install, configure, and verify HOL Guard runtime protection for supported AI agent hosts, including tool execution, secret access, prompt injection, and extension supply-chain checks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# HOL Guard Agent Runtime Security

## Overview

Deploy and verify HOL Guard as a local runtime policy layer for supported agent tools, file access, package operations, and MCP integrations. Establish which events the selected host actually intercepts, configure scoped policy, and collect evidence from controlled allowed and denied actions. The upstream “antivirus” description is product positioning, not a guarantee against compromise or prompt injection.

## When to Use This Skill

- Install or diagnose HOL Guard for a supported agent host.
- Configure restrictions on sensitive files, destructive tools, package installation, or changed agent extensions.
- Verify whether an MCP call or shell operation is blocked before execution.
- Investigate unexpected approvals, missing receipts, false positives, or coverage changes after an upgrade.
- Review plugins and skills while distinguishing static findings from runtime enforcement.

## Core Capabilities

1. **Establish the enforcement boundary.** Record the HOL Guard release, host/version, operating system, architecture, workspace, and requested surfaces. Read the version-matched installation guide and support matrix; record the documentation revision when available. Map each surface to native hooks, managed MCP proxy, launch inspection, prompt screening, or post-tool observation. Treat documentation as declared capability until a live host test demonstrates enforcement.

   Examples from the upstream support matrix, to recheck for the installed release:

   | Host | Declared boundary | Important limitation |
   | --- | --- | --- |
   | Codex | Native `PreToolUse` and `PermissionRequest`; managed launch checks hooks | Hook presence alone does not prove live blocking |
   | Claude Code | Native pre-tool and permission hooks | Guard does not install `UserPromptSubmit` interception |
   | Cursor | Shell, MCP, pre-tool, and file-read hooks | Native approval prompts alone are not Guard evidence |
   | OpenCode | Guard pre-tool plugin and managed MCP path | Coverage depends on the plugin and managed routes being active |
   | Cline | Native hooks or AgentPlugin transport | Native `PostToolUse` is observation-only; JetBrains needs live pre-tool proof |

   Other hosts require their own matrix entries. Do not infer support from a repository topic or from another adapter's behavior. The matrix documents fail-open hook crash/timeout behavior for hosts including Kimi, Grok, and ZCode; a fail-closed internal runtime does not eliminate host-level fail-open behavior.

2. **Install and integrate within the requested scope.** Verify Python 3.10+, pipx, and a compatible native runtime against the release installation guide. Preserve existing host configuration and a recoverable pre-change copy before installation. Select an approved release rather than silently upgrading an existing deployment.

   Basic upstream installation and diagnostic commands:

   ```bash
   pipx install hol-guard
   hol-guard init
   hol-guard --version
   hol-guard status
   ```

   For a Codex deployment, the documented explicit path is:

   ```bash
   hol-guard install codex
   hol-guard doctor codex
   hol-guard run codex --dry-run
   hol-guard run codex
   ```

   Use the selected adapter's documented identifier and commands for another host. The dry run records artifact state; it is not a live tool-denial test. Keep optional cloud connection outside the deployment unless requested; local protection does not require an account. Do not install software merely to answer a review-only request.

3. **Configure and validate policy.** Inspect effective policy and its precedence before changing it. Express the requested allow, review, or block behavior using the installed release's supported schema and matchers. Scope exceptions to the needed operation or artifact; preserve required protections. GuardPolicy `v1alpha1` is an alpha interchange contract, not a stable-schema promise.

   For releases supporting the documented policy interface:

   ```bash
   hol-guard policy export --output ./policy-before.yaml --json
   hol-guard policy validate ./guard-policy.yaml --json
   hol-guard policy fmt ./guard-policy.yaml --check --json
   hol-guard policy diff ./guard-policy.yaml --json
   hol-guard policy import ./guard-policy.yaml --merge --dry-run --json
   hol-guard policy check ./guard-policy.yaml --command "docker ps"
   ```

   Inspect the candidate diff and predicted decision before applying an authorized local change with `hol-guard policy import ./guard-policy.yaml --merge --apply --json`. Export and validate effective policy afterward; compare it to the intended merged result, accounting for retained rules. Downloading a recipe or validating a YAML file does not activate it. For managed signed bundles, follow the operator runbook's capability, signature, acknowledgement, and rollback gates; do not enable canonical enforcement simply to make a smoke test pass.

4. **Review supply-chain surfaces.** Use `hol-guard inventory`, `hol-guard diff codex` for the Codex example, and `hol-guard supply-chain scan` to inspect tracked changes and dependency findings. Use `hol-guard supply-chain sync` when advisory refresh is in scope, and record advisory freshness and ecosystem gaps. Confirm which package-manager operations and MCP routes the selected integration intercepts.

   If static plugin inspection is requested, the separate `plugin-scanner` distribution supports `plugin-scanner scan ./candidate --format json`. A scan result, quality grade, clean advisory lookup, or inventoried MCP server does not establish runtime protection or absence of malicious behavior. Do not execute an untrusted plugin or install a suspicious package to validate a finding.

5. **Prove allowed and denied behavior safely.** Use an isolated disposable workspace with synthetic content and no production credentials. Define expected decisions from the active policy before testing. First inspect command classification without execution:

   ```bash
   hol-guard command test 'rm -rf ./build'
   hol-guard command explain 'git clean -ndx'
   ```

   These commands classify text without executing the embedded command or creating an approval. They do not prove the host consumes a denial. Then run a small live test matrix through the selected host:

   | Case | Controlled fixture | Evidence required |
   | --- | --- | --- |
   | Allowed tool | Read a benign fixture or perform an allowed harmless command | Successful result plus correlated Guard decision |
   | Denied tool/file access | A harmless action explicitly denied by the test policy, or a protected file containing only dummy text | Host denial, Guard receipt, and confirmation the action did not occur or content was withheld |
   | Approval path, if required | Harmless action configured for review | Pending request and no execution before the authorized decision |
   | Prompt/result screening, if supported | Synthetic instruction attempting access to a dummy protected fixture | Actual screen/rewrite/block evidence; identify when observation occurs after execution |
   | MCP, if in scope | Harmless call on a trusted test server through the managed route | Matching tool identity, decision, and suppressed call for the deny case |

   Review `hol-guard receipts` and `hol-guard approvals`; redact sensitive values. A pending approval is not a final denial. For adapters that encode denial in JSON while exiting zero, inspect the envelope and host behavior rather than treating exit status as the verdict. Do not auto-approve a test request to clear a queue. Mark unsupported or untested cases explicitly.

6. **Report limits and recover from failures.** Distinguish a declared capability, installed configuration, synthetic classifier test, and live enforcement proof. Post-tool observation cannot undo an executed action. Direct tool/server access outside managed routes, disabled integrations, omitted host events, host error handling, and unrecognized command or package forms can leave gaps; maintain host permissions and sandbox controls independently.

   Stop rollout on malformed or unsupported policy, signature/hash failure, unexpected allow, missing live evidence, or an unexplained policy mismatch. Preserve diagnostic evidence, restore the prior known-good policy/configuration using the installed release's rollback procedure, and repeat the benign allow/deny checks. Do not weaken policy or disable hooks to obtain a passing result. Assess false positives with command explanations and narrowly scoped policy changes.

## Inputs / Outputs

**Inputs:** Target host and version; OS/architecture; workspace and configuration scope; approved HOL Guard release; desired policy decisions; MCP/package/plugin inventory; and disposable synthetic fixtures. Existing configuration, effective policy, and the installed release's support matrix are the baseline.

**Outputs:** Installation and integration status; reviewed policy diff and effective-policy identity; a per-surface coverage table; redacted allowed/denied test evidence with receipt references; supply-chain findings where requested; and recovery instructions with known gaps.

**Completion criterion:** Within the 15-minute target, one selected host has a verified integration, the intended effective policy, and at least one live allowed action and one live denied action with evidence of non-execution or withholding. Report additional surfaces as verified, unsupported, failed, or not tested. If prerequisites or live testing prevent completion, deliver a bounded incomplete report with the exact blocker; do not label installation or a static scan as successful runtime verification.

## References

- [HOL Guard source repository](https://github.com/hashgraph-online/hol-guard)
- [Installation and platform support](https://github.com/hashgraph-online/hol-guard/blob/main/docs/guard/get-started.md)
- [Host support matrix, failure behavior, and explicit exclusions](https://github.com/hashgraph-online/hol-guard/blob/main/docs/guard/harness-support.md)
- [GuardPolicy format and validation commands](https://github.com/hashgraph-online/hol-guard/blob/main/spec/guard-policy/v1alpha1/README.md)
- [Policy validation, import, rollout, and stop conditions](https://github.com/hashgraph-online/hol-guard/blob/main/spec/guard-policy/v1alpha1/operator-runbook.md)
- [Policy recipes and activation boundary](https://github.com/hashgraph-online/hol-guard/blob/main/docs/guard/policy-recipes.md)
- [Policy rollback runbook](https://github.com/hashgraph-online/hol-guard/blob/main/spec/guard-policy/v1alpha1/rollback-runbook.md)
