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

# Add the 'Skills' directory to path so we can import our new agent skills
# In a real package, this would be handled by pip install
sys.path.append(os.path.abspath("../../../"))

from Agentic_AI.Multi_Agent_Systems.debate_supervisor import Supervisor

class MedicalSupervisor(Supervisor):
    """
    A specialized debate supervisor for Clinical Diagnosis.
    Inherits from the generic Multi-Agent Supervisor.
    """
    def __init__(self):
        super().__init__()
        # Override agents with medical personas
        self.agents["Optimist"].persona = "You are an aggressive Oncologist. You prioritize catching every potential cancer early, even if it means aggressive testing."
        self.agents["Critic"].persona = "You are a conservative Palliative Care specialist. You prioritize patient quality of life and avoiding over-treatment/unnecessary biopsies."

    def diagnose(self, patient_case):
        print(f"--- Medical Board: Discussing Case ---")
        print(f"Case: {patient_case}")
        summary = self.conduct_debate(patient_case, rounds=2)
        print("\n--- Final Recommendation ---")
        print(summary)

if __name__ == "__main__":
    case = "75M with 2cm lung nodule, stable for 2 years, but now has mild cough. Heavy smoker history. COPD."
    
    med_board = MedicalSupervisor()
    med_board.diagnose(case)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
