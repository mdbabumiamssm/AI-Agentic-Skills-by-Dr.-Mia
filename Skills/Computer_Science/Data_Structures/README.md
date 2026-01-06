# Vector Database (In-Memory)

**Category:** Computer Science / Data Structures
**Difficulty:** Intermediate

## Overview
A Vector Database is a specialized database optimized for storing and querying high-dimensional vectors. It is the backbone of modern Semantic Search and RAG applications.

While standard databases (SQL) query for exact matches (`WHERE id = 5`), Vector Databases query for **similarity** (`WHERE distance(vector, query) is minimized`).

## Implementation Details (`vector_store.py`)
This is a "Flat Index" implementation.
- **Storage:** Python Dictionary (Hash Map) for O(1) retrieval by ID.
- **Search:** Brute-force Euclidean Distance calculation across all $N$ vectors.
- **Complexity:** $O(N \cdot D)$ per query.

## Production Alternatives
In real-world systems (like Faiss or Annoy), brute force is too slow for millions of vectors. They use Approximate Nearest Neighbor (ANN) algorithms:
1.  **HNSW (Hierarchical Navigable Small World):** Graph-based.
2.  **IVF (Inverted File Index):** Clustering-based.
3.  **PQ (Product Quantization):** Compression-based.
