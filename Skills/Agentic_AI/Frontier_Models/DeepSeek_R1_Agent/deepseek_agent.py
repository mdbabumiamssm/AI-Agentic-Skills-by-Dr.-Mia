import time
from typing import Dict, Any

class DeepSeekAgent:
    """
    Wrapper for DeepSeek R1, focusing on exposing the Chain-of-Thought.
    """
    def __init__(self):
        self.model = "deepseek-r1"

    def solve(self, problem: str) -> Dict[str, Any]:
        print(f"🐳 [DeepSeek R1] Pondering: {problem[:50]}...")
        
        # Simulated Chain of Thought
        cot = (
            "<think>\n"
            "1. User wants a proof for infinitely many primes of form 4n+3.\n"
            "2. Recall Euclid's theorem for primes.\n"
            "3. Assume finite set {p1, ..., pk}.\n"
            "4. Construct N = 4*p1*...*pk - 1.\n"
            "5. Analyze factors of N modulo 4.\n"
            "6. Contradiction found.\n"
            "</think>"
        )
        
        answer = "The proof follows by contradiction. Consider N = 4(p1...pk) - 1. N is of the form 4n+3..."
        
        return {
            "model": self.model,
            "reasoning": cot,
            "answer": answer
        }

if __name__ == "__main__":
    agent = DeepSeekAgent()
    res = agent.solve("Solve integral of x^2 * e^x dx")
    print(res["reasoning"])

