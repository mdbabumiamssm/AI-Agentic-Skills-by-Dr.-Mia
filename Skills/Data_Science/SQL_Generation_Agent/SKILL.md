# SQL Generation Agent Skill

**Domain:** Data Science / Database Management
**Status:** Active

## Overview
The SQL Generation Agent provides a structured interface for translating natural language into SQL queries. This is a common "text-to-SQL" pattern utilized in modern data platforms. It supports providing schema context (DDL) to ground the language model and improve query accuracy.

## Capabilities
- `generate_sql(natural_language_query, schema_context)`: Generates a SQL query string based on the user's intent and an optional database schema.

## Implementation Details
- **Language:** Python
- **Dependencies:** Standard library
- **Main Class:** `SQLGenerationAgent`
- **Note:** The current implementation uses a heuristic mock for demonstration. In production, this integrates with the project's core LLM adapters.
