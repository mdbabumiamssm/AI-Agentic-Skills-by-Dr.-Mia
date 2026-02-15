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

# React & Next.js Performance Rules

## 1. Eliminate Waterfalls
*   **Rule**: Avoid sequential data fetching in parent-child relationships where the child waits for the parent to finish fetching before it starts.
*   **Fix**: Use `Promise.all` for parallel fetching at the route level (Page/Layout), or use the Preload pattern.
*   **Context**: In Next.js App Router, prefer fetching data in Server Components.

## 2. Optimize Bundle Size
*   **Rule**: Do not import heavy libraries (e.g., large charting libs, 3D renderers) in the initial bundle if they are not immediately visible.
*   **Fix**: Use `next/dynamic` (Next.js) or `React.lazy` (React) with `suspense` to lazy load these components.
*   **Example**:
    ```tsx
    import dynamic from 'next/dynamic'
    const HeavyChart = dynamic(() => import('./HeavyChart'), { loading: () => <p>Loading...</p> })
    ```

## 3. Server vs. Client Components
*   **Rule**: Default to Server Components. Only use `"use client"` when you need:
    *   Interactivity (onClick, onChange).
    *   React State (useState, useReducer).
    *   Lifecycle effects (useEffect).
    *   Browser-only APIs (window, localStorage).
*   **Benefit**: Reduces the amount of JavaScript sent to the client.

## 4. Image Optimization
*   **Rule**: Avoid standard `<img>` tags for local assets or known remote sources.
*   **Fix**: Use `next/image`.
    *   Always define `width` and `height` (or `fill`).
    *   Use `sizes` prop for responsive images to serve the correct variant.

## 5. Prevent Unnecessary Re-renders
*   **Rule**: Don't pass inline objects or arrow functions as props to memoized child components unless necessary.
*   **Fix**: Use `useMemo` for objects/arrays and `useCallback` for functions *only if* the child component is wrapped in `React.memo` or is expensive to render. Premature optimization here adds overhead.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->