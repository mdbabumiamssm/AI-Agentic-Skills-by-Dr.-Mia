import os
import sys
import asyncio
import json
from typing import Dict, Any, List, Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

# --- Import Resolution ---
current_dir = os.path.dirname(os.path.abspath(__file__))
# Root project dir
project_root = os.path.abspath(os.path.join(current_dir, "../../"))
if project_root not in sys.path:
    sys.path.append(project_root)

# Platform dir
platform_dir = os.path.join(project_root, "platform")
if platform_dir not in sys.path:
    sys.path.append(platform_dir)

# Core Systems
try:
    from adapters.runtime_adapter import llm
    from Skills.Research_Tools.Literature_Mining.mining_agent import LiteratureMiningAgent
    from Skills.Drug_Discovery.Molecule_Design.evolution_agent import MoleculeEvolutionAgent
    from Skills.Clinical.Safety.safety_agent import SafetyOfficerAgent
    from Skills.Clinical.Trial_Matching.trial_matching_agent import TrialMatchingAgent
except ImportError as e:
    print(f"Critical Import Error: {e}")
    # We allow the server to start even if agents fail, but endpoints will error
    LiteratureMiningAgent = None
    MoleculeEvolutionAgent = None
    SafetyOfficerAgent = None
    TrialMatchingAgent = None

app = FastAPI(title="BioKernel Enterprise", version="2026.3.0-PRO")

class AgentRequest(BaseModel):
    query: str
    mode: str = "direct" # 'direct' or 'workflow'

class PatientProfile(BaseModel):
    id: str
    age: int
    diagnosis: str
    status: str
    ecog_score: int

class BioKernel:
    def __init__(self):
        self.miner = LiteratureMiningAgent() if LiteratureMiningAgent else None
        self.chemist = MoleculeEvolutionAgent() if MoleculeEvolutionAgent else None
        self.deputy = SafetyOfficerAgent() if SafetyOfficerAgent else None
        self.matcher = TrialMatchingAgent() if TrialMatchingAgent else None

    async def execute_discovery_workflow(self, initial_query: str):
        """
        Orchestrates: Mining -> Evolution -> Safety
        """
        if not (self.miner and self.chemist and self.deputy):
            return "Error: One or more agents unavailable."

        print(f"\n🚀 [Kernel] Initializing Discovery Workflow: '{initial_query}'")
        results = []

        # 1. Mine
        mining_result = self.miner.mine_for_targets(initial_query)
        target = mining_result.get("target_protein", "Unknown")
        results.append(f"1. Mining: Found target '{target}' ({mining_result.get('status')})")

        if mining_result.get("status") == "found":
            # 2. Design
            drug_result = self.chemist.evolve(target)
            candidate = drug_result.get("top_candidate")
            results.append(f"2. Design: Evolved candidate '{candidate}'")

            # 3. Safety
            safety_result = self.deputy.inspect_output(f"Proposed drug {candidate}")
            results.append(f"3. Safety: Status is {safety_result.get('status').upper()}")
        
        return "\n".join(results)

    async def execute_trial_match(self, patient: Dict[str, Any]):
        if not self.matcher:
            return "Error: Trial Matcher unavailable."
        
        matches = self.matcher.match_patient(patient)
        return {"patient_id": patient.get("id"), "matches": matches}

kernel = BioKernel()

@app.get("/")
async def health_check():
    return {"status": "active", "system": "BioKernel OS"}

@app.post("/v1/workflow/discovery")
async def run_discovery(request: AgentRequest):
    report = await kernel.execute_discovery_workflow(request.query)
    return {"report": report}

@app.post("/v1/workflow/trial_match")
async def run_trial_match(profile: PatientProfile):
    result = await kernel.execute_trial_match(profile.dict())
    return result

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)