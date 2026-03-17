# SQLite MCP Agent Skill

**Domain:** Agentic AI / Model Context Protocol (MCP)
**Status:** Active
**Trend Context:** 2026 MCP Standardization

## Overview
The SQLite MCP (Model Context Protocol) Agent provides a standardized interface for Large Language Models to connect to, introspect, and safely query local SQLite databases. As the "Agent Internet" protocols like MCP gain massive traction in 2026, skills like this are fundamental for giving models secure access to local structured data.

## Capabilities
- `get_schema()`: Introspects the database and returns a JSON representation of tables and columns to build prompt context.
- `execute_query(query)`: Safely executes `SELECT` statements and formats the result set for LLM consumption. Enforces read-only operations for security.

## Implementation Details
- **Language:** Python
- **Dependencies:** Standard library (`sqlite3`)
- **Main Class:** `SQLiteMCPAgent`
