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
name: 'core-python-best-practices'
description: 'Essential guidelines for writing modern, type-safe, and idiomatic Python 3 code.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
  - write_file
---


# Core Python Best Practices

This skill defines the coding standards for Python development within the project. It emphasizes modern features, type safety, and readability.

## When to Use This Skill

*   **New Scripts**: Starting a new agent or tool.
*   **Refactoring**: Modernizing legacy code.
*   **Library Design**: Creating reusable modules.

## Core Capabilities

1.  **Type Hinting**: Mandatory use of `typing` module or native types (Python 3.9+).
2.  **Data Classes**: Using `@dataclass` or `Pydantic` for data containers instead of raw dictionaries/tuples.
3.  **Modern Control Flow**: Using `match/case` (Python 3.10) where appropriate.
4.  **Error Handling**: Proper use of `try/except` chains and custom exceptions.

## Workflow

1.  **Define Interface**: Start with function signatures and type hints.
2.  **Select Structure**: Choose between a simple function, a class, or a dataclass.
3.  **Implement**: Write logic using list comprehensions and generators where possible.
4.  **Document**: Add docstrings (Google or NumPy style).

## Example Usage

**User**: "Write a function to process a list of users."

**Agent Action**:
1.  Reads `references/rules.md`.
2.  Generates:
    ```python
    from dataclasses import dataclass
    
    @dataclass
    class User:
        id: int
        name: str

    def process_users(users: list[User]) -> None:
        """Processes a list of users."""
        for user in users:
            ...
    ```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->