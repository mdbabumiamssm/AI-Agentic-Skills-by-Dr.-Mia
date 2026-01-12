import os
import json
import asyncio
import time
from typing import Dict, Any, List, Optional
from fastapi import FastAPI
from pydantic import BaseModel

# Import Core Systems
try:
    from Skills.Computer_Science.Distributed_Systems.event_bus import bus
    from Skills.Agentic_AI.Memory_Systems.blackboard import shared_context
    # Dynamic Imports for Workflow Agents
    from Skills.Research_Tools.Literature_Mining.mining_agent import LiteratureMiningAgent
    from Skills.Drug_Discovery.Molecule_Design.evolution_agent import MoleculeEvolutionAgent
    from Skills.Clinical.Safety.safety_agent import SafetyOfficerAgent
except ImportError:
    bus = None
    shared_context = None

app = FastAPI(title="BioKernel Enterprise", version="2026.3.0-PRO")

class AgentRequest(BaseModel):
    query: str
    mode: str = "direct" # 'direct' or 'workflow'

class BioKernel:
    def __init__(self):
        self.miner = LiteratureMiningAgent()
        self.chemist = MoleculeEvolutionAgent()
        self.deputy = SafetyOfficerAgent()

    async def execute_workflow(self, initial_query: str):
        """
        Enterprise Workflow Orchestration:
        1. Mining: Literature search for targets.
        2. Design: Generative chemistry for candidates.
        3. Compliance: Safety verification.
        """
        print(f"\n🚀 [Kernel] Initializing Workflow: '{initial_query}'")
        results = []

        # --- Step 1: Literature Mining ---
        if bus: await bus.publish("step_start", {"step": "Mining"}, "Kernel")
        mining_result = self.miner.mine_for_targets(initial_query)
        shared_context.write("latest_discovery", mining_result, "Miner")
        results.append(f"📰 Miner found: {mining_result.get('target_protein', 'Nothing')}")
        await asyncio.sleep(1) 

        if mining_result.get("status") == "found":
            target = mining_result["target_protein"]
            
            # --- Step 2: Evolution (Triggered by Shared Context) ---
            if bus: await bus.publish("step_start", {"step": "Evolution"}, "Kernel")
            print(f"  ↪️ Triggering Designer for target: {target}")
            drug_result = self.chemist.evolve(target)
            candidate = drug_result["top_candidate"]
            shared_context.write("candidate_drug", candidate, "Designer")
            results.append(f"🧬 Designed candidate: {candidate} (Score: {drug_result['score']:.2f})")
            await asyncio.sleep(1)

            # --- Step 3: Safety Check (Triggered by Shared Context) ---
            if bus: await bus.publish("step_start", {"step": "Safety"}, "Kernel")
            print(f"  ↪️ Triggering Safety Officer to inspect: {candidate}")
            safety_result = self.deputy.inspect_output(f"Proposed drug {candidate} for target {target}")
            shared_context.write("safety_report", safety_result, "Safety")
            
            status_icon = "✅" if safety_result["status"] == "approved" else "❌"
            results.append(f"{status_icon} Safety check: {safety_result['status']}")

        return "\n".join(results)

kernel = BioKernel()

@app.post("/v1/workflow/run")
async def run_workflow(request: AgentRequest):
    report = await kernel.execute_workflow(request.query)
    return {"report": report}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001) # Port 8001 for Workflow Engine
