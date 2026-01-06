"Agentic Memory System (2026 Update)

Integrates:
1. Short-term Episodic Memory (Rolling window)
2. Long-term Semantic Memory (Vector Database)
3. Real Embeddings (SentenceTransformer)
"

import sys
import os
import time
import json
import numpy as np

# Add path to Computer Science utils
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../Computer_Science/Data_Structures')))

try:
    from vector_store import VectorStore
    HAS_VECTOR_STORE = True
except ImportError:
    HAS_VECTOR_STORE = False

# Try import embeddings
try:
    from sentence_transformers import SentenceTransformer
    EMBEDDING_MODEL = SentenceTransformer('all-MiniLM-L6-v2')
    DIMENSION = 384
    HAS_EMBEDDINGS = True
except ImportError:
    HAS_EMBEDDINGS = False
    DIMENSION = 5 # Mock dimension

class MemorySystem:
    def __init__(self, agent_name="Agent"):
        self.agent_name = agent_name
        self.episodic = [] # List of dicts
        
        if HAS_VECTOR_STORE:
            self.semantic = VectorStore(dimension=DIMENSION)
        else:
            self.semantic = None
            print("Warning: VectorStore dependency missing.")

    def _get_embedding(self, text):
        if HAS_EMBEDDINGS:
            return EMBEDDING_MODEL.encode(text).tolist()
        else:
            # Deterministic mock embedding for consistent testing without GPU/Libs
            # seeded by length and char codes
            np.random.seed(len(text)) 
            return np.random.rand(DIMENSION).tolist()

    def log_interaction(self, role: str, content: str):
        """Logs an event to Episodic Memory."""
        event = {
            "timestamp": time.time(),
            "role": role,
            "content": content
        }
        self.episodic.append(event)
        
        # Heuristic: If it's a user statement, save to Semantic Memory as a 'Fact'
        # In a real system, an LLM would extract 'claims' first.
        if role == "user" and self.semantic:
            doc_id = f"mem_{len(self.episodic)}"
            vec = self._get_embedding(content)
            self.semantic.add(doc_id, vec, content, metadata={"source": "conversation"})

    def retrieve_context(self, query: str, k: int = 2) -> str:
        """
        Retrieves a hybrid context:
        1. Recent messages (Episodic)
        2. Relevant past facts (Semantic)
        """
        # 1. Episodic (Last 3 turns)
        recent_history = "\n".join(
            [f"{e['role']}: {e['content']}" for e in self.episodic[-3:]]
        )
        
        # 2. Semantic
        semantic_context = ""
        if self.semantic:
            query_vec = self._get_embedding(query)
            results = self.semantic.search(query_vec, k=k)
            if results:
                facts = [f"- {r[0]['text']}" for r in results]
                semantic_context = "Relevant Knowledge:\n" + "\n".join(facts)
        
        return f"{semantic_context}\n\nRecent Chat:\n{recent_history}"

if __name__ == "__main__":
    mem = MemorySystem()
    print("--- Initializing Memory System ---")
    
    # 1. User teaches the agent
    print("User: My favorite color is blue.")
    mem.log_interaction("user", "My favorite color is blue.")
    mem.log_interaction("agent", "I've noted that.")
    
    # 2. Distractor
    print("User: Let's talk about math.")
    mem.log_interaction("user", "What is 2+2?")
    mem.log_interaction("agent", "It is 4.")
    
    # 3. Recall Test
    query = "What does the user like?"
    print(f"\n--- RECALL TEST: '{query}' ---")
    context = mem.retrieve_context(query)
    print(context)