# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""
Foundational Vector Store (2026 Update)

A lightweight, persistent, in-memory vector database optimized for:
1. RAG Systems
2. Semantic Memory
3. Educational demonstration of vector search algorithms

Features:
- Cosine Similarity (standard for semantic search)
- Metadata Filtering
- Persistence (JSON/Numpy serialization)
"""

import numpy as np
import json
import os
import heapq
from typing import List, Dict, Tuple, Optional

class VectorDocument:
    def __init__(self, id: str, vector: np.ndarray, text: str, metadata: Dict):
        self.id = id
        self.vector = vector
        self.text = text
        self.metadata = metadata

class VectorStore:
    def __init__(self, dimension: int, persist_dir: Optional[str] = None):
        self.dimension = dimension
        self.documents: Dict[str, VectorDocument] = {}
        self.persist_dir = persist_dir
        
        if persist_dir and os.path.exists(os.path.join(persist_dir, "index.json")):
            self.load()

    def add(self, id: str, vector: List[float], text: str, metadata: Dict = None):
        if len(vector) != self.dimension:
            raise ValueError(f"Vector dim {len(vector)} != {self.dimension}")
        
        # Normalize for Cosine Similarity optimization (L2 Norm)
        vec_np = np.array(vector, dtype=np.float32)
        norm = np.linalg.norm(vec_np)
        if norm > 0:
            vec_np = vec_np / norm
            
        self.documents[id] = VectorDocument(id, vec_np, text, metadata or {})

    def search(self, query_vector: List[float], k: int = 3, filter_dict: Dict = None) -> List[Tuple[Dict, float]]:
        """
        Returns top-k documents sorted by Cosine Similarity.
        Score range: [-1, 1] (1 is identical).
        """
        query_np = np.array(query_vector, dtype=np.float32)
        norm = np.linalg.norm(query_np)
        if norm > 0:
            query_np = query_np / norm
            
        scores = []
        
        for doc_id, doc in self.documents.items():
            # 1. Metadata Filtering
            if filter_dict:
                match = True
                for key, val in filter_dict.items():
                    if doc.metadata.get(key) != val:
                        match = False
                        break
                if not match:
                    continue
            
            # 2. Cosine Similarity (Dot product of normalized vectors)
            score = float(np.dot(query_np, doc.vector))
            scores.append((doc, score))
        
        # Top-K
        # We use nlargest because higher cosine similarity is better
        top_k = heapq.nlargest(k, scores, key=lambda x: x[1])
        
        return [
            ({
                "id": doc.id,
                "text": doc.text,
                "metadata": doc.metadata
            }, score) 
            for doc, score in top_k
        ]

    def save(self):
        if not self.persist_dir:
            return
        
        os.makedirs(self.persist_dir, exist_ok=True)
        
        # Save Metadata & Text
        index_data = {
            doc_id: {
                "text": doc.text, 
                "metadata": doc.metadata,
                "vector": doc.vector.tolist()
            } 
            for doc_id, doc in self.documents.items()
        }
        
        with open(os.path.join(self.persist_dir, "index.json"), "w") as f:
            json.dump(index_data, f)
            
    def load(self):
        path = os.path.join(self.persist_dir, "index.json")
        with open(path, "r") as f:
            data = json.load(f)
            
        for doc_id, val in data.items():
            self.add(doc_id, val["vector"], val["text"], val["metadata"])

if __name__ == "__main__":
    # Demo
    vs = VectorStore(dimension=3)
    
    # "King" - "Man" + "Woman" approx "Queen"
    vs.add("1", [1.0, 0.0, 0.0], "King", {"category": "royal"})
    vs.add("2", [0.9, 0.1, 0.0], "Queen", {"category": "royal"})
    vs.add("3", [0.0, 1.0, 0.0], "Apple", {"category": "fruit"})
    
    print("Searching for vector close to King...")
    results = vs.search([0.95, 0.05, 0.0], k=2)
    
    for doc, score in results:
        print(f"[{score:.4f}] {doc['text']} ({doc['metadata']})")
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
