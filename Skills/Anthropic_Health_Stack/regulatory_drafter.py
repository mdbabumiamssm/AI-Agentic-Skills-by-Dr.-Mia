import time
import sys
import os
import json
from typing import Dict, Any, Optional, List

# Adjust path to find platform module if running standalone
if __name__ == "__main__":
    # Add project root
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
    # Add platform directory to avoid 'platform' stdlib collision
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../platform")))

try:
    from optimizer.meta_prompter import PromptOptimizer, ModelTarget
except ImportError:
    # Fallback if paths are tricky
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../platform")))
    from optimizer.meta_prompter import PromptOptimizer, ModelTarget

# Import SelfCorrectionAgent
try:
    from Skills.Agentic_AI.Agent_Architectures.Self_Correction.self_correction_agent import SelfCorrectionAgent
except ImportError:
    # Fallback for relative imports
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../Agentic_AI/Agent_Architectures/Self_Correction")))
    from self_correction_agent import SelfCorrectionAgent

class RegulatoryDrafter:
    """
    Regulatory Agent: Regulatory Drafter
    
    Demonstrates:
    1. Integration with Meta-Prompter for optimized prompt engineering.
    2. Agentic Self-Correction Loop for high-quality drafting.
    3. Chain-of-Thought transparency (<thinking> blocks).
    4. Structured Output.
    """
    
    def __init__(self, model: str = "claude-3-5-sonnet-20241022", adapter=None):
        self.model = model
        self.optimizer = PromptOptimizer()
        self.adapter = adapter
        self.corrector = SelfCorrectionAgent()

    def set_adapter(self, adapter):
        """Inject a real LLM adapter (e.g. from platform.adapters)."""
        self.adapter = adapter

    def draft_submission(self, 
                         section_name: str,
                         clinical_data: str, 
                         guidance_text: str,
                         style_guide: Optional[str] = None) -> Dict[str, Any]:
        """
        Drafts a regulatory submission section using an agentic feedback loop.
        
        Args:
            section_name: e.g. "Pediatric Assessment Waiver", "Nonclinical Overview"
            clinical_data: The raw data/study results.
            guidance_text: Relevant FDA/EMA guidance.
            style_guide: Optional formatting/tone instructions.
            
        Returns:
            Dictionary with status, trace, and draft content.
        """
        
        # 1. Define the Task for the Agent
        task = f"""
        Role: Senior Regulatory Affairs Writer
        Task: Draft the '{section_name}' section for a submission.
        
        Context:
        - Clinical Data: {clinical_data}
        - Guidance: {guidance_text}
        """
        
        if style_guide:
            task += f"\n- Style Guide: {style_guide}"
            
        task += """
        
        Requirements:
        - Analyze the data against the guidance.
        - Cite specific regulations/guidance sections where applicable.
        - Maintain a formal, objective regulatory tone.
        - Output format: Markdown.
        """
        
        # 2. Define Success Criteria for the Critic (Self-Correction)
        criteria = [
            "Strictly adhere to the provided guidance text.",
            "Use formal, objective regulatory language.",
            "Explicitly reference data points from the clinical data.",
            "Ensure the structure matches standard CTD formats.",
            "No hallucinations; only use provided data."
        ]

        print(f"--- RegulatoryDrafter: Processing '{section_name}' with Agentic Loop ---")
        
        # 3. Execution (Real or Mock Agent Loop)
        # The SelfCorrectionAgent handles the loop. 
        # In a real scenario, SelfCorrectionAgent would call an LLM. 
        # Here, SelfCorrectionAgent runs a simulation if no LLM is wired, 
        # or we can rely on its internal mock if adapter is missing.
        
        result = self.corrector.run_cycle(task, criteria, max_iterations=2)

        return {
            "status": "success",
            "model": self.model,
            "prompt_used": task, # The initial task prompt
            "draft": result["final_output"],
            "iterations": result["iterations"],
            "history": result["history"]
        }

if __name__ == "__main__":
    agent = RegulatoryDrafter()
    
    # Test 1: Pediatric Waiver
    print(">> TEST 1: Pediatric Waiver")
    res1 = agent.draft_submission(
        section_name="Pediatric Assessment Waiver",
        clinical_data="Indication: Wet AMD. Mean age 72. No pediatric incidence.",
        guidance_text="FDA Guidance: Waivers granted if disease does not occur in children."
    )
    print("\n[FINAL DRAFT]:")
    print(res1["draft"])
    
    # Test 2: Nonclinical Overview
    print("\n>> TEST 2: Nonclinical Overview")
    res2 = agent.draft_submission(
        section_name="Nonclinical Overview",
        clinical_data="Ki = 0.5 nM. No off-target effects in Cerep panel.",
        guidance_text="ICH M3(R2): Nonclinical safety studies for conduct of human clinical trials."
    )
    print("\n[FINAL DRAFT]:")
    print(res2["draft"])