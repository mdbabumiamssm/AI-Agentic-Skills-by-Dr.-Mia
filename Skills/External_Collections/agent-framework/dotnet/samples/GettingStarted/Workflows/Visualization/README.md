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

﻿# Workflow Visualization Sample

This sample demonstrates how to visualize workflows using `ToMermaidString()` and `ToDotString()` extension methods. It uses a map-reduce workflow with fan-out/fan-in patterns as an example.

## Running the Sample

```bash
dotnet run
```

## Output Formats

The sample generates two visualization formats:

### Mermaid
Paste the output into any Mermaid-compatible viewer (GitHub, Mermaid Live Editor, etc.):

![Mermaid Visualization](Resources/mermaid_render.png)

### DOT (Graphviz)
Render with Graphviz (requires `graphviz` to be installed):

```bash
dotnet run | tail -n +20 | dot -Tpng -o workflow.png
```

![Graphviz Visualization](Resources/graphviz_render.png)

## Usage

```csharp
Workflow workflow = BuildWorkflow();

// Generate Mermaid format
string mermaid = workflow.ToMermaidString();

// Generate DOT format
string dotString = workflow.ToDotString();
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->