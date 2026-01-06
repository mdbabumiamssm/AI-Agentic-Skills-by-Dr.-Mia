"""
MedPrompt Utilities (2026)

This module implements advanced prompt engineering strategies specifically designed
for clinical and biomedical tasks, inspired by "MedPrompt" and "Chain-of-Verification".

Key Strategies:
1. Dynamic Few-Shot Selection (kNN based - mocked here)
2. Chain-of-Thought (CoT) with Ensemble Refinement
3. Structuring outputs to FHIR standards
"""

from typing import List, Dict, Optional

class MedPromptEngine:
    def __init__(self, model_name: str = "gpt-4-turbo-biomed"):
        self.model_name = model_name

    def generate_chain_of_thought_prompt(self, task: str, clinical_text: str) -> str:
        """
        Constructs a prompt that forces the model to reason step-by-step 
        before providing the final answer.
        """
        prompt = f"""
You are an expert clinical AI assistant.
Task: {task}

Input Clinical Text:
"""
{clinical_text}
"""

Instructions:
1. Analyze the patient's demographics and chief complaint.
2. Identify positive and negative findings in the text.
3. Reason through the differential diagnosis or treatment plan.
4. ONLY after reasoning, provide the final structured output.

Output Format:
[Reasoning]
... your step-by-step logic ...

[Final Output]
... the requested summary/extraction ...
"""
        return prompt

    def chain_of_verification(self, initial_response: str, clinical_text: str) -> str:
        """
        Generates a verification prompt to error-check the initial response.
        """
        prompt = f"""
We have an initial analysis of a clinical note. Your job is to VERIFY it for hallucinations or omissions.

Original Note:
"""
{clinical_text}
"""

Initial Analysis:
"""
{initial_response}
"""

Verification Steps:
1. List all medical facts in the Initial Analysis.
2. Check each fact against the Original Note.
3. If a fact is NOT in the note, mark it as HALLUCINATION.
4. If a critical finding (e.g., severe vitals) is missing, mark as OMISSION.

Return the corrected analysis.
"""
        return prompt

    def generate_ensemble_prompts(self, task: str, clinical_text: str, num_paths: int = 5) -> List[str]:
        """
        Generates multiple distinct prompts to encourage diverse reasoning paths 
        (Self-Consistency).
        """
        prompts = []
        strategies = [
            "Think step-by-step focusing on ruling out dangerous conditions first.",
            "Think step-by-step focusing on chronological events.",
            "Think step-by-step focusing on evidence-based guidelines.",
            "Act as a skeptic and challenge every assumption.",
            "Act as a concise summarizer focusing on billing codes."
        ]
        
        for i in range(num_paths):
            strategy = strategies[i % len(strategies)]
            prompt = f"""
Task: {task}
Strategy: {strategy}
Input: {clinical_text}
Output: Provide your best answer.
"""
            prompts.append(prompt)
        return prompts

    def consensus_vote(self, responses: List[str]) -> str:
        """
        Placeholder for a majority-voting or aggregator logic.
        In production, this would use embeddings to find the centroid answer.
        """
        # Simple length-based heuristic for now (preferring more detailed answers)
        return max(responses, key=len)


    def format_as_fhir_json(self, summary_data: Dict) -> str:
        """
        Helper to wrap data in a pseudo-FHIR Bundle structure.
        """
        import json
        
        fhir_bundle = {
            "resourceType": "Bundle",
            "type": "document",
            "entry": []
        }
        
        # Example mapping
        if "patient_age" in summary_data:
            fhir_bundle["entry"].append({
                "resource": {
                    "resourceType": "Patient",
                    "id": "generated-id",
                    "extension": [{
                        "url": "http://hl7.org/fhir/StructureDefinition/patient-age",
                        "valueInteger": summary_data["patient_age"]
                    }]
                }
            })
            
        if "diagnosis" in summary_data:
             fhir_bundle["entry"].append({
                "resource": {
                    "resourceType": "Condition",
                    "code": {
                        "text": summary_data["diagnosis"]
                    }
                }
            })
            
        return json.dumps(fhir_bundle, indent=2)

if __name__ == "__main__":
    # Demo Usage
    engine = MedPromptEngine()
    
    note = "Pt 45M c/o chest pain. ECG shows ST elevation. Troponin +."
    
    print("--- 1. Chain of Thought Prompt ---")
    print(engine.generate_chain_of_thought_prompt("Summarize", note))
    
    print("\n--- 2. Chain of Verification Prompt ---")
    fake_response = "Patient is 45F with abdominal pain." # Intentional error
    print(engine.chain_of_verification(fake_response, note))
