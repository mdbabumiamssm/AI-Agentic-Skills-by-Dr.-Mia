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
name: 'moai-spec-first-agentic-devkit'
description: 'Use MOAI ADK SPEC-first agentic workflows to turn requirements into specs, TDD checks, DDD models, agent-team plans, and multilingual delivery.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# MOAI Spec-First Agentic Development Kit

## Overview

This skill guides agents through MOAI ADK-style SPEC-first development: convert user intent into explicit requirements, tests, domain models, implementation plans, and quality gates before coding. Use it to operationalize agentic software delivery where requirements clarity, TDD, DDD, and coordinated AI agent roles matter more than ad hoc code generation.

MOAI ADK is a Go-based CLI and workflow package for Claude Code that presents agent teams, skills, scaffolding, and multilingual delivery patterns for spec-driven projects. This skill adapts those practices into a concise agent workflow without assuming a specific repository layout unless the target project already uses MOAI ADK.

## When to Use This Skill

- The user asks for SPEC-first, spec-driven, requirements-first, or EARS-style software development.
- The project needs TDD or quality gates before implementation.
- The task benefits from DDD concepts such as bounded contexts, aggregates, entities, value objects, or domain services.
- The user wants AI agent team planning, role assignment, or multi-agent orchestration for a software project.
- A repository uses or evaluates `modu-ai/moai-adk`, `moai`, Claude Code agents, or MOAI ADK skills.
- The deliverable spans multiple programming or documentation languages and needs a consistent specification-to-implementation workflow.

## Core Capabilities

1. **Requirement capture** - Translate user goals into clear behavioral requirements, acceptance criteria, constraints, and open questions before proposing implementation.
2. **SPEC-first planning** - Structure work around a specification artifact that defines scope, expected behavior, non-goals, interfaces, and quality gates.
3. **EARS-style acceptance criteria** - Express critical behavior using event-driven, state-driven, optional, unwanted, and ubiquitous requirement patterns when useful.
4. **TDD workflow design** - Identify tests that should exist before code changes, including unit, integration, contract, regression, and smoke tests.
5. **DDD modeling** - Map the domain into bounded contexts, ubiquitous language, aggregates, entities, value objects, repositories, and services when the problem has domain complexity.
6. **Agent team orchestration** - Break large work into bounded agent roles such as product/spec, domain modeling, test design, implementation, review, and documentation.
7. **Quality gate execution** - Define and run checks for formatting, linting, tests, type checks, security-sensitive paths, docs, and release readiness.
8. **Multilingual delivery alignment** - Keep requirements, tests, implementation, and documentation synchronized across languages or generated project scaffolds.

## Inputs / Outputs

**Inputs**

- User goal, feature request, bug report, or product requirement.
- Existing repository files, issue text, design notes, or architecture documentation.
- Target language, framework, runtime, deployment environment, and test command constraints.
- Any existing MOAI ADK configuration, generated agents, skills, specs, or scaffolding.

**Outputs**

- A concise specification with requirements, assumptions, non-goals, acceptance criteria, and quality gates.
- A test-first implementation plan that identifies the minimum useful tests before code changes.
- A DDD model when domain structure is relevant.
- An agent-role work breakdown for larger tasks.
- Concrete commands or file edits needed to implement, verify, and document the requested change.
- A final verification summary showing which quality gates passed, failed, or were not run.

## References

- MOAI ADK GitHub repository: https://github.com/modu-ai/moai-adk
