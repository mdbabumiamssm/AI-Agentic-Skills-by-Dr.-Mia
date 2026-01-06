# Dijkstra's Algorithm

**Category:** Computer Science / Algorithms / Graphs
**Difficulty:** Intermediate

## Overview
Dijkstra's algorithm finds the shortest paths between nodes in a graph, which may represent, for example, road networks. It was conceived by computer scientist Edsger W. Dijkstra in 1956 and published three years later.

It is a **Greedy Algorithm** that uses a **Priority Queue** to explore the most promising (shortest path) nodes first.

## Complexity
- **Time Complexity:** $O((V + E) \log V)$ using a binary heap (priority queue).
- **Space Complexity:** $O(V)$ to store distances and predecessors.

Where $V$ is the number of vertices and $E$ is the number of edges.

## Implementation Details
We use Python's built-in `heapq` module to implement the priority queue efficiently.

### Key Concepts
1.  **Relaxation:** The process of checking if a path through a neighbor is shorter than the currently known path.
    ```python
    if dist[u] + weight(u, v) < dist[v]:
        dist[v] = dist[u] + weight(u, v)
    ```
2.  **Priority Queue:** Ensures we always process the unvisited node with the smallest known distance.

## Usage
Run the standalone script to see the algorithm in action on a sample graph.

```bash
python dijkstra.py
```

## Applications
- Google Maps (Routing)
- Network Routing Protocols (OSPF)
- IP Routing
