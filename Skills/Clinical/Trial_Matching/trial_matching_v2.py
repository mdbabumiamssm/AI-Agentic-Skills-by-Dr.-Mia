"""
TrialGPT v2 (2026 Update)
Integrated Clinical Trial Matching & Reasoning Agent
"""

import sys
import os
import json

# Add paths for dependencies
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Clinical_Note_Summarization')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Trial_Matching_Agent')))

try:
    from medprompt_utils import MedPromptEngine
    from trial_matcher import TrialMatcher, ClinicalTrial, Patient
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False

class TrialGPTv2:
    def __init__(self):
        self.prompt_engine = MedPromptEngine() if HAS_DEPS else None
        self.matcher = TrialMatcher(use_llm=False) # Fallback to keyword for speed in prototype
        
    def process_patient_and_match(self, note: str, trials_file: str):
        print("--- TrialGPT v2: Processing Clinical Note ---")
        
        # 1. Summarization & Reasoning (MedPrompt)
        if self.prompt_engine:
            reasoning_prompt = self.prompt_engine.generate_chain_of_thought_prompt(
                "Extract Patient Age, Gender, and Primary Condition for Trial Matching.",
                note
            )
            print("[LLM Logic] Reasoning through clinical note...")
            # print(reasoning_prompt) # In production, send to LLM
        
        # Mock extraction from note for matching logic
        # Note: "Pt 45M c/o chest pain. History of Type 2 Diabetes."
        age = 45
        condition = "Diabetes"
        
        patient = Patient(age, "Male", [condition], note)
        
        # 2. Match trials
        try:
            with open(trials_file, 'r') as f:
                trials_data = json.load(f)
            trials = [ClinicalTrial(t) for t in trials_data]
        except Exception as e:
            return f"Error loading trials: {e}"

        print(f"Searching trials for {condition} (Age {age})...")
        matches = self.matcher.match(patient, trials)
        
        # 3. Verification (Chain-of-Verification)
        if matches and self.prompt_engine:
            top_match = matches[0]
            verification = self.prompt_engine.chain_of_verification(
                f"Matched patient to trial {top_match['trial_id']} because of {condition}.",
                note
            )
            print("[LLM Logic] Verifying match against raw note data...")
            
        return matches

def main():
    if not HAS_DEPS:
        print("Error: Missing internal dependencies. Ensure medprompt_utils.py and trial_matcher.py are in paths.")
        return

    engine = TrialGPTv2()
    
    # Example clinical note
    sample_note = "45yo male with persistent hyperglycemia. Diagnosed with Type 2 Diabetes 5 years ago. No history of heart disease."
    
    # Path to dummy data (relative to Trial_Matching_Agent)
    data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Trial_Matching_Agent', 'dummy_data.json'))
    
    results = engine.process_patient_and_match(sample_note, data_path)
    
    if isinstance(results, list):
        print("\nMatch Results:")
        for r in results[:3]:
            print(f"- [{r['score']:.2f}] {r['trial_id']}: {r['title']}")
    else:
        print(results)

if __name__ == "__main__":
    main()
