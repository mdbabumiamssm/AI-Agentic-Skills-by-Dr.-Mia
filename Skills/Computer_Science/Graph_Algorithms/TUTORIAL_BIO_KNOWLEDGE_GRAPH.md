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

# Tutorial: Building a Biomedical Knowledge Graph

**Domain:** Computer Science / Graph Theory
**Focus:** GraphRAG & Drug Repurposing
**Level:** Intermediate

## Introduction

Biology is a network, not a list.
*   Drugs bind to Proteins.
*   Proteins regulate Pathways.
*   Pathways cause Diseases.

If you just search text (RAG), you miss these connections. **GraphRAG** combines text search with graph traversal to answer questions like "How does Drug X affect Disease Y?" even if no single document mentions both.

## The Data Structure

In `knowledge_graph.py`, we implement a simple graph:
*   **Nodes:** Entities (e.g., "Aspirin", "COX-2").
*   **Edges:** Relationships (e.g., "inhibits").

## Algorithm: Shortest Path (Dijkstra)

Why do we need this?
Imagine you want to know *why* a drug works.
The algorithm finds the chain:
`Drug A` -> inhibits -> `Protein B` -> regulates -> `Process C` -> causes -> `Disease D`.

```python
kg = BioKnowledgeGraph()
# ... add nodes and edges ...
path = kg.find_shortest_path("Drug A", "Disease D")
# Output: ['Drug A', 'Protein B', 'Process C', 'Disease D']
```

## Algorithm: Subgraph Context (for LLMs)

When a user asks "Tell me about Gene X", passing the whole database to the LLM is too big.
We use `get_subgraph_context("Gene X", depth=1)` to grab *only* the immediate neighbors.

```python
context = kg.get_subgraph_context("BCR-ABL")
# Returns: ["(Imatinib) --[inhibits]--> (BCR-ABL)", "(BCR-ABL) --[causes]--> (CML)"]
```
You feed these triples to the LLM as "Context".

## Assignments

1.  **Build:** Add a new node "Gleevec" (synonym for Imatinib) and link it.
2.  **Analyze:** Use `calculate_centrality()` to find the most important protein in your graph (the one with the most connections).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->