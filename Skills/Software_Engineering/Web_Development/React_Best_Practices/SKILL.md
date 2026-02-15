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

---
name: 'react-best-practices'
description: 'A comprehensive guide and rule set for writing clean, performant, and maintainable React code.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
  - write_file
---


# React Best Practices

This skill provides a set of architectural and coding standards for React applications. It is designed to be used by agents when generating, refactoring, or reviewing React code to ensure alignment with modern industry standards.

## When to Use This Skill

*   **Code Generation**: When asked to create new React components or hooks.
*   **Code Review**: When analyzing existing React code for anti-patterns.
*   **Refactoring**: When optimizing performance or improving readability.
*   **Migration**: When upgrading from class components to functional components + hooks.

## Core Capabilities

1.  **Component Architecture**: Enforces separation of concerns (Container/Presentational) or Composition patterns.
2.  **Hooks Usage**: Validates rules of hooks, custom hook extraction, and dependency array correctness.
3.  **State Management**: Guidance on `useState` vs `useReducer` vs Context vs External Stores (Zustand/Redux).
4.  **Performance Optimization**: Identifies unnecessary re-renders and suggests `useMemo`/`useCallback` usage appropriately.

## Workflow

1.  **Analyze Context**: Determine if the task involves a single component or a complex feature.
2.  **Consult Rules**: Refer to `references/rules.md` for specific guidelines on naming, file structure, and typing.
3.  **Generate/Refactor**: Apply the rules to the code.
4.  **Validate**: Ensure no "prop drilling" or "huge components" are introduced.

## Example Usage

**User**: "Create a user profile component that fetches data."

**Agent Action**:
1.  Reads `references/rules.md` to check data fetching patterns (e.g., "Use React Query or SWR over useEffect for data").
2.  Generates the component:
    ```tsx
    // user-profile.tsx
    import { useQuery } from '@tanstack/react-query';
    
    // ... implementation following rules ...
    ```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->