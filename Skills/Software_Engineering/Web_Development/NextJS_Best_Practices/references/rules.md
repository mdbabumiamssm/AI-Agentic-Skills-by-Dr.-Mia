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

# Next.js (App Router) Best Practices Rules

## 1. Project Structure (App Router)

*   **Colocation**: Keep components used in a single route inside that route's folder (unless shared).
*   **Private Folders**: Use `_folderName` to exclude folders from routing.
*   **Route Groups**: Use `(groupName)` to organize routes without affecting the URL structure.

## 2. Server vs. Client Components

*   **Default to Server**: Server Components reduce bundle size and access backend resources directly.
*   **"use client" Boundary**: Place the directive at the top of the file only when you need:
    *   `useState`, `useEffect`
    *   Event listeners (`onClick`, `onChange`)
    *   Browser-only APIs (`window`, `localStorage`)
*   **Leaf Pattern**: Push client components down the tree. Don't make the root layout a client component.

## 3. Data Fetching & Caching

*   **Native Fetch**: Use the extended `fetch` API.
    ```typescript
    // Revalidates every hour
    fetch('https://...', { next: { revalidate: 3600 } })
    ```
*   **No Store**: For dynamic data, use `{ cache: 'no-store' }`.
*   **Parallel Data Fetching**: Initiate requests in parallel to prevent waterfalls.
    ```typescript
    const [user, posts] = await Promise.all([getUser(), getPosts()]);
    ```

## 4. Server Actions

*   **Forms**: Use `action` prop on forms to invoke server-side logic directly.
*   **Validation**: Use **Zod** to validate inputs inside the server action.
*   **Revalidation**: Use `revalidatePath` or `revalidateTag` to update cached data after mutation.

## 5. SEO & Metadata

*   **Static Metadata**: Export `metadata` object from `layout.tsx` or `page.tsx`.
*   **Dynamic Metadata**: Export `generateMetadata({ params })` function for dynamic routes.
*   **Open Graph**: Use `opengraph-image.tsx` for dynamic OG image generation.

## 6. Image Optimization

*   **next/image**: Always use the `<Image>` component instead of `<img>`.
*   **Sizing**: Provide `width` and `height` or `fill` to prevent Layout Shift (CLS).
*   **Priority**: Add `priority` prop to the LCP (Largest Contentful Paint) image (usually the hero image).

## 7. Security

*   **Taint**: Use `experimental_taintUniqueValue` to prevent sensitive data from being passed to Client Components.
*   **Server-Only**: Install `server-only` package to ensure sensitive modules are not imported into client code.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->