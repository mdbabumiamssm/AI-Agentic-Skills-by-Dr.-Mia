<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Core Python Best Practices Rules

## 1. Type Hinting (Typing)

*   **Strict Typing**: All function arguments and return values MUST have type hints.
*   **Native Types**: Use `list[str]`, `dict[str, int]`, `tuple[int, ...]` (Python 3.9+) instead of `List`, `Dict`, `Tuple` imports from `typing` where possible.
*   **Optional**: Use `Optional[str]` or `str | None` (Python 3.10+) for nullable values.
*   **Any**: Avoid `Any`. If strictly necessary, explain why in a comment.

## 2. Data Structures

*   **Dataclasses**: Use `@dataclass` for classes that primarily store data.
    ```python
    from dataclasses import dataclass

    @dataclass
    class Point:
        x: float
        y: float
    ```
*   **Pydantic**: Use `Pydantic` models when data validation (e.g., API payloads, JSON parsing) is required.

## 3. Documentation

*   **Docstrings**: All public functions and classes must have a docstring.
*   **Style**: Prefer Google Style or NumPy Style.
    ```python
    def calculate_distance(p1: Point, p2: Point) -> float:
        """Calculates Euclidean distance between two points.

        Args:
            p1: The first point.
            p2: The second point.

        Returns:
            The distance as a float.
        """
    ```

## 4. Modern Idioms

*   **f-strings**: Always use f-strings `f"{val}"` instead of `%` formatting or `.format()`.
*   **Pathlib**: Use `pathlib.Path` instead of `os.path` for file system manipulation.
*   **Context Managers**: Use `with` statements for file I/O, locks, and database connections.

## 5. Error Handling

*   **Specific Exceptions**: Catch specific exceptions (`ValueError`, `KeyError`) rather than bare `except:`.
*   **Custom Exceptions**: Define custom exception classes for domain-specific errors.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->