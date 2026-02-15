# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import heapq

def dijkstra(graph, start_node):
    """
    Implements Dijkstra's algorithm to find the shortest path from start_node to all other nodes.
    
    Args:
        graph: Dict[str, Dict[str, float]] 
               Adjacency list where graph[u][v] = weight
        start_node: str
    
    Returns:
        distances: Dict[str, float] (shortest distance to each node)
        path: Dict[str, str] (predecessor map to reconstruct paths)
    """
    # Priority queue stores (distance, node)
    pq = [(0, start_node)]
    
    # Distances initialized to infinity
    distances = {node: float('inf') for node in graph}
    distances[start_node] = 0
    
    # To reconstruct the path
    predecessors = {node: None for node in graph}
    
    while pq:
        current_dist, current_node = heapq.heappop(pq)
        
        # If we found a shorter path previously, skip
        if current_dist > distances[current_node]:
            continue
        
        for neighbor, weight in graph[current_node].items():
            distance = current_dist + weight
            
            # Only consider this new path if it's better
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                predecessors[neighbor] = current_node
                heapq.heappush(pq, (distance, neighbor))
                
    return distances, predecessors

def get_shortest_path(predecessors, target_node):
    """Reconstructs the path from start to target."""
    path = []
    current = target_node
    while current is not None:
        path.append(current)
        current = predecessors[current]
    return path[::-1]  # Reverse to get start -> target

if __name__ == "__main__":
    # Example Graph
    # A --1--> B --2--> C
    # |        |
    # 4        1
    # v        v
    # D --1--> E
    
    graph = {
        'A': {'B': 1, 'D': 4},
        'B': {'C': 2, 'E': 1},
        'C': {},
        'D': {'E': 1},
        'E': {}
    }
    
    start = 'A'
    dists, preds = dijkstra(graph, start)
    
    print(f"Shortest distances from {start}: {dists}")
    print(f"Path to E: {get_shortest_path(preds, 'E')}")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
