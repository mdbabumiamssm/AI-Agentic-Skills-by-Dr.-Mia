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
name: sql-generation-agent
description: Generates and evaluates grounded Text-to-SQL queries using schema context, production monitoring, and constrained decoding patterns.
keywords:
  - sql-generation
  - text-to-sql
  - data-science
  - schema-grounding
  - constrained-decoding
measurable_outcome: Generate schema-grounded SQL and evaluate production Text-to-SQL outputs with documented semantic checks.
license: Proprietary
compatibility:
  - system: Python 3.10+
allowed-tools:
  - read_file
  - run_shell_command
---

# SQL Generation Agent Skill

**Domain:** Data Science / Database Management
**Status:** Active

## Overview
The SQL Generation Agent provides a structured interface for translating natural language into SQL queries. This is a common "text-to-SQL" pattern utilized in modern data platforms. It supports providing schema context (DDL) to ground the language model and improve query accuracy.

## Capabilities
- `generate_sql(natural_language_query, schema_context)`: Generates a SQL query string based on the user's intent and an optional database schema.
- `evaluate_sql_production(user_question, enriched_reformulation, generated_sql)`: Applies STEF-style production Text-to-SQL evaluation by extracting semantic specifications from the natural-language inputs and generated SQL, aligning normalized features, and producing schema-agnostic scoring and monitoring without reference SQL.
- `template_constrained_decoding(recurring_workload, user_question, schema_context)`: For recurring Text-to-SQL workloads, mine reusable templates from historical NL-SQL pairs, reject unmatched queries, select matched templates with natural language inference, and enforce grammar-constrained SQL generation for safer low-latency production behavior.

## Implementation Details
- **Language:** Python
- **Dependencies:** Standard library
- **Main Class:** `SQLGenerationAgent`
- **Note:** The current implementation uses a heuristic mock for demonstration. In production, this integrates with the project's core LLM adapters.

## References
- http://arxiv.org/abs/2604.28049v1
- http://arxiv.org/abs/2604.28028v1
