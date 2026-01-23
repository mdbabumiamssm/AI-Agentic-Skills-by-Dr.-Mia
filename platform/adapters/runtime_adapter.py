from typing import Dict, Any, List, Optional
import time

class RuntimeLLMAdapter:
    """
    Unified Runtime Adapter for executing LLM prompts.
    
    Modes:
    - 'mock': Returns heuristic/template-based responses (default).
    - 'api': (Future) Connects to actual Claude/OpenAI APIs.
    """
    
    def __init__(self, mode: str = "mock", api_key: Optional[str] = None):
        self.mode = mode
        self.api_key = api_key

    def complete(self, system: str, user_prompt: str, model: str = "default") -> str:
        """
        Generates a text completion.
        """
        if self.mode == "mock":
            return self._mock_response(system, user_prompt)
        else:
            raise NotImplementedError("API mode not yet implemented.")

    def _mock_response(self, system: str, user_prompt: str) -> str:
        """
        Generates context-aware mock responses to make demos feel alive.
        """
        p_lower = user_prompt.lower()
        
        # --- Literature Mining Scenarios ---
        if "target" in p_lower and "mine" in p_lower:
            return "Analysis of recent literature (Nature, Cell, 2024-2025) reveals a novel target 'GPRC5D' implicated in resistant multiple myeloma. Confidence: High."
        
        # --- Chemistry / Evolution Scenarios ---
        if "smiles" in p_lower or "evolve" in p_lower:
            return "Proposed modification: Add a fluorine group to the ortho-position of the phenyl ring to improve metabolic stability. New SMILES: CC(=O)OC1=C(F)C=CC=C1C(=O)O"
            
        # --- Safety Scenarios ---
        if "inspect" in p_lower or "safety" in p_lower:
            if "ricin" in p_lower or "poison" in p_lower:
                return "CRITICAL ALERT: Dangerous content detected. The input references known toxins. Request REJECTED."
            return "Safety Review: PASSED. No harmful compounds, PII, or hallucinated medical claims detected."

        # --- Clinical Trial Scenarios ---
        if "trial" in p_lower or "match" in p_lower:
            return "Patient matches Trial NCT04561234 (Phase 3). Inclusion criteria met: Age > 18, Histology confirmed. Exclusion criteria cleared."

        # Default
        return f"[MockLLM] Processed query: {user_prompt[:50]}..."

# Singleton instance for easy import
llm = RuntimeLLMAdapter()
