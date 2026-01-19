import sys
import os
import json
from typing import Dict, Any, Optional, List

# Adjust path to find platform module
if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../platform")))

try:
    from optimizer.meta_prompter import PromptOptimizer, ModelTarget
except ImportError:
    # Fallback for when running from different contexts
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../platform")))
    from optimizer.meta_prompter import PromptOptimizer, ModelTarget

class SelfCorrectionAgent:
    """
    Implements a 'Reflexion' or 'Self-Correction' pattern.
    
    Workflow:
    1. Generator: Produces an initial draft based on a task.
    2. Critic: Reviews the draft against specific criteria or ground truth.
    3. Refiner: Updates the draft based on the critic's feedback.
    """

    def __init__(self, model_target: ModelTarget = ModelTarget.OPENAI):
        self.optimizer = PromptOptimizer()
        self.model_target = model_target

    def run_cycle(self, task: str, criteria: List[str], max_iterations: int = 2) -> Dict[str, Any]:
        """
        Executes the self-correction cycle.
        """
        history = []
        
        # Step 1: Initial Draft
        draft_prompt = self._build_generator_prompt(task)
        # In a real system, this would call the LLM. We will simulate/mock for the skill demo.
        current_draft = self._mock_llm_call(draft_prompt, "draft")
        
        history.append({"step": "initial_draft", "content": current_draft})
        print(f"[{self.model_target.name}] generated initial draft.")

        for i in range(max_iterations):
            # Step 2: Critique
            critique_prompt = self._build_critic_prompt(task, current_draft, criteria)
            critique = self._mock_llm_call(critique_prompt, "critique")
            history.append({"step": f"critique_{i+1}", "content": critique})
            
            if "NO_ISSUES" in critique:
                print(f"Critique passed at iteration {i+1}.")
                break
                
            # Step 3: Refine
            refine_prompt = self._build_refiner_prompt(task, current_draft, critique)
            current_draft = self._mock_llm_call(refine_prompt, "refine")
            history.append({"step": f"refined_draft_{i+1}", "content": current_draft})
            print(f"Refined draft at iteration {i+1}.")

        return {
            "final_output": current_draft,
            "history": history,
            "iterations": i + 1
        }

    def _build_generator_prompt(self, task: str) -> str:
        raw = f"Task: {task}\nGenerate a comprehensive solution."
        return self.optimizer.optimize(raw, self.model_target)

    def _build_critic_prompt(self, task: str, draft: str, criteria: List[str]) -> str:
        criteria_str = "\n- ".join(criteria)
        return f"""
        Task: {task}
        Current Draft: {draft}
        
        Critique this draft based on:
        - {criteria_str}
        
        If it meets all criteria, output 'NO_ISSUES'. Otherwise, list specific actionable improvements.
        """

    def _build_refiner_prompt(self, task: str, draft: str, critique: str) -> str:
        return f"""
        Task: {task}
        Draft: {draft}
        Critique: {critique}
        
        Rewrite the draft to address the critique.
        """

    def _mock_llm_call(self, prompt: str, mode: str) -> str:
        """
        Simulates LLM responses for demonstration purposes.
        In a production BioKernel env, this would use `platform.adapters`.
        """
        prompt_lower = prompt.lower()
        
        # Domain: Prior Authorization / Appeals
        if "appeal" in prompt_lower or "denial" in prompt_lower:
            if mode == "draft":
                return "To Whom It May Concern, I am writing to appeal the denial for Patient X. They have diabetes."
            elif mode == "critique":
                if "contraindication" not in prompt_lower and "policy" not in prompt_lower:
                    return "Critique: The draft is too generic. It fails to mention the specific contraindication (renal failure) or the payer policy code."
                return "NO_ISSUES"
            elif mode == "refine":
                return "Dear Medical Reviewer, This is an appeal for Patient X. The denial states 'step therapy required'. However, the patient has Stage 4 CKD (eGFR 25), which is a contraindication for Metformin per Policy 101. Therefore, the requested medication is medically necessary."

        # Domain: Regulatory
        elif "regulatory" in prompt_lower or "waiver" in prompt_lower:
            if mode == "draft":
                return "We request a waiver. Kids don't get this disease."
            elif mode == "critique":
                if "cfr" not in prompt_lower:
                    return "Critique: Too informal. Needs to cite 21 CFR 314.55."
                return "NO_ISSUES"
            elif mode == "refine":
                return "Pursuant to 21 CFR 314.55(c)(2), we request a full waiver. The disease etiology in pediatric populations is nonexistent."

        # Default / Medical Triage
        else:
            if mode == "draft":
                return "The patient has a fever. We should give Tylenol."
            elif mode == "critique":
                if "dosage" not in prompt_lower and "history" not in prompt_lower: 
                    return "The draft fails to mention dosage and checking patient history for liver disease."
                return "NO_ISSUES"
            elif mode == "refine":
                return "The patient has a fever. Check liver history. If clear, administer 650mg Acetaminophen q6h."
        
        return " [ LLM Output ] "

if __name__ == "__main__":
    agent = SelfCorrectionAgent()
    task = "Recommend treatment for 38C fever in adult."
    criteria = ["Include Dosage", "Check Contraindications (Liver)", "Specify Frequency"]
    
    result = agent.run_cycle(task, criteria)
    print(json.dumps(result, indent=2))
