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

# Everything Server - Extension Points

**[Architecture](architecture.md)
| [Project Structure](structure.md)
| [Startup Process](startup.md)
| [Server Features](features.md)
| Extension Points
| [How It Works](how-it-works.md)**

## Adding Tools

- Create a new file under `tools/` with your `registerXTool(server)` function that registers the tool via `server.registerTool(...)`.
- Export and call it from `tools/index.ts` inside `registerTools(server)`.

## Adding Prompts

- Create a new file under `prompts/` with your `registerXPrompt(server)` function that registers the prompt via `server.registerPrompt(...)`.
- Export and call it from `prompts/index.ts` inside `registerPrompts(server)`.

## Adding Resources

- Create a new file under `resources/` with your `registerXResources(server)` function using `server.registerResource(...)` (optionally with `ResourceTemplate`).
- Export and call it from `resources/index.ts` inside `registerResources(server)`.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->