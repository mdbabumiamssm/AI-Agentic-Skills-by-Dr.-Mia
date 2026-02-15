# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import sys
import os
import json
from typing import Dict, Any, List

# Adjust path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../../.."))
platform_dir = os.path.join(project_root, "platform")

if platform_dir not in sys.path:
    sys.path.append(platform_dir)

from adapters.runtime_adapter import llm

class TrialMatchingAgent:
    """
    Agent for matching patients to clinical trials using intelligent criteria mapping.
    """
    
    def __init__(self):
        # Mock Database of Clinical Trials
        self.trial_db = [
            {
                "nct_id": "NCT04561234",
                "title": "Phase 3 Study of Novel CAR-T in Multiple Myeloma",
                "phase": "Phase 3",
                "criteria": "Inclusion: Age > 18, Relapsed/Refractory MM, ECOG 0-1. Exclusion: Active infection, CNS involvement."
            },
            {
                "nct_id": "NCT09876543",
                "title": "Study of XYZ-123 in Solid Tumors",
                "phase": "Phase 1",
                "criteria": "Inclusion: Advanced solid tumors, progressive disease. Exclusion: Prior immunotherapy."
            }
        ]

    def match_patient(self, patient_profile: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Finds matching trials for a given patient profile.
        """
        print(f"🏥 [Matcher] Matching patient {patient_profile.get('id', 'Unknown')}...")
        
        matches = []
        patient_str = json.dumps(patient_profile)
        
        for trial in self.trial_db:
            # Ask LLM to evaluate eligibility
            prompt = f"""
            Patient: {patient_str}
            Trial: {json.dumps(trial)}
            
            Task: Determine if the patient is eligible for this trial.
            Output: "ELIGIBLE" or "NOT ELIGIBLE" followed by a brief reason.
            """
            
            decision = llm.complete("You are a clinical research coordinator.", prompt)
            
            if "ELIGIBLE" in decision or "matches" in decision.lower():
                matches.append({
                    "trial_id": trial["nct_id"],
                    "title": trial["title"],
                    "reasoning": decision
                })
                
        return matches

if __name__ == "__main__":
    agent = TrialMatchingAgent()
    patient = {
        "id": "PT-001",
        "age": 65,
        "diagnosis": "Multiple Myeloma",
        "status": "Relapsed",
        "ecog_score": 1
    }
    results = agent.match_patient(patient)
    print(json.dumps(results, indent=2))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
