import time
import sys
import os
from typing import Dict, Any, Optional

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

class RegulatoryDrafter:
    """
    Anthropic Health Stack Agent: Regulatory Drafter
    
    Demonstrates:
    1. Integration with Meta-Prompter for optimized prompt engineering.
    2. Long Context Handling (up to 200k tokens).
    3. Chain-of-Thought transparency (<thinking> blocks).
    4. Structured Output.
    """
    
    def __init__(self, model: str = "claude-3-5-sonnet-20241022", adapter=None):
        self.model = model
        self.optimizer = PromptOptimizer()
        self.adapter = adapter

    def set_adapter(self, adapter):
        """Inject a real LLM adapter (e.g. from platform.adapters)."""
        self.adapter = adapter

    def draft_submission(self, 
                         section_name: str,
                         clinical_data: str, 
                         guidance_text: str,
                         style_guide: Optional[str] = None) -> Dict[str, Any]:
        """
        Drafts a regulatory submission section based on data and guidance.
        
        Args:
            section_name: e.g. "Pediatric Assessment Waiver", "Nonclinical Overview"
            clinical_data: The raw data/study results.
            guidance_text: Relevant FDA/EMA guidance.
            style_guide: Optional formatting/tone instructions.
            
        Returns:
            Dictionary with status, trace, and draft content.
        """
        
        # 1. Optimize Prompt
        base_context = f"""
        Role: Senior Regulatory Affairs Writer
        Task: Draft the '{section_name}' section for a submission.
        
        Context:
        - Clinical Data: {clinical_data}
        - Guidance: {guidance_text}
        """
        
        if style_guide:
            base_context += f"\n- Style Guide: {style_guide}"
            
        raw_prompt = f"""
        {base_context}
        
        Requirements:
        - Analyze the data against the guidance.
        - Cite specific regulations/guidance sections where applicable.
        - Maintain a formal, objective regulatory tone.
        - Output format: Markdown.
        """
        
        optimized_prompt = self.optimizer.optimize(raw_prompt, ModelTarget.CLAUDE)
        
        print(f"--- RegulatoryDrafter: Processing '{section_name}' with {self.model} ---")
        
        # 2. Execution (Real or Mock)
        if self.adapter:
            # In a real BioKernel setup, we'd use the adapter
            response = self.adapter.generate(optimized_prompt)
            return {
                "status": "success",
                "model": self.model,
                "prompt_used": optimized_prompt,
                "draft": response
            }
        else:
            return self._run_simulation(section_name, optimized_prompt)

    def _run_simulation(self, section_name: str, prompt: str) -> Dict[str, Any]:
        """Mock simulation for demonstration/testing without API keys."""
        print(">> Claude is thinking...")
        time.sleep(0.8)
        
        # Dynamic-ish thinking block
        thinking_process = f"""
<thinking>
1.  **Analyze the Request**: Draft '{section_name}'.
2.  **Review Guidance**: Checking provided guidance for constraints on '{section_name}'.
3.  **Review Data**: Synthesizing clinical data points relevant to this section.
4.  **Formulate Argument**: Constructing narrative flow: Introduction -> Data Presentation -> Conclusion.
5.  **Drafting Strategy**: Adhering to CTD format.
</thinking>
"""
        
        # Placeholder drafts for common demos
        if "Pediatric" in section_name:
            draft_text = """
**Draft Submission: Pediatric Assessment Waiver Request**

Pursuant to 21 CFR 314.55(c)(2), the Sponsor requests a full waiver...

**Rationale:**
The indication is adult-onset (e.g., AMD) and does not exist in pediatric populations.
"""
        elif "Nonclinical" in section_name:
             draft_text = """
**2.4 Nonclinical Overview**

**2.4.1 Introduction**
The nonclinical program for [Drug X] was designed to support chronic administration...

**2.4.2 Pharmacology**
Primary pharmacology studies demonstrated high affinity binding (Ki = 0.5 nM)...
"""
        else:
            draft_text = f"**{section_name}**\n\n[Draft content generated based on input data...]"

        return {
            "status": "success",
            "model": self.model,
            "prompt_used": prompt,
            "trace": thinking_process.strip(),
            "draft": draft_text.strip()
        }

if __name__ == "__main__":
    agent = RegulatoryDrafter()
    
    # Test 1: Pediatric Waiver
    print(">> TEST 1: Pediatric Waiver")
    res1 = agent.draft_submission(
        section_name="Pediatric Assessment Waiver",
        clinical_data="Indication: Wet AMD. Mean age 72.",
        guidance_text="FDA Guidance on Pediatric Studies..."
    )
    print(res1["draft"][:100] + "...")
    
    # Test 2: Nonclinical Overview
    print("\n>> TEST 2: Nonclinical Overview")
    res2 = agent.draft_submission(
        section_name="Nonclinical Overview",
        clinical_data="Ki = 0.5 nM. No off-target effects.",
        guidance_text="ICH M3(R2)"
    )
    print(res2["draft"][:100] + "...")