import sys
import os
from typing import Dict, Any

# Adjust path to find platform module
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../../../.."))
platform_dir = os.path.join(project_root, "platform")

if platform_dir not in sys.path:
    sys.path.append(platform_dir)

from adapters.runtime_adapter import llm

class LiteratureMiningAgent:
    """
    Agent responsible for mining scientific literature for novel targets.
    Uses LLM (simulated or real) to extract insights.
    """
    def __init__(self):
        pass

    def mine_for_targets(self, query: str) -> Dict[str, Any]:
        print(f"📰 [Miner] Scanning literature for '{query}'...")
        
        # Ask the LLM to perform the mining
        prompt = f"Mine recent literature for targets related to: {query}"
        insight = llm.complete("You are a biomedical researcher.", prompt)
        
        # Parse the 'insight' (in a real app, we'd request JSON)
        # Here we mock the parsing based on the mock response
        
        if "novel target" in insight or "found" in insight:
            # simple extraction heuristic for demo
            return {
                "status": "found",
                "target_protein": "GPRC5D" if "GPRC5D" in insight else "Unknown-Target",
                "summary": insight,
                "confidence": 0.95
            }
        
        return {"status": "none", "reason": "No high-confidence targets found."}

if __name__ == "__main__":
    agent = LiteratureMiningAgent()
    print(agent.mine_for_targets("multiple myeloma"))
