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
    from Skills.Agentic_AI.Memory_Systems.blackboard import town_square
    # Dynamic Imports for "The Posse"
    from Skills.Research_Tools.Literature_Mining.mining_agent import LiteratureMiningAgent
    from Skills.Drug_Discovery.Molecule_Design.evolution_agent import MoleculeEvolutionAgent
    from Skills.Clinical.Safety.safety_agent import SafetyOfficerAgent
except ImportError:
    bus = None
    town_square = None

app = FastAPI(title="BioKernel Frontier", version="2026.3.0-POSSE")

class AgentRequest(BaseModel):
    query: str
    mode: str = "direct" # 'direct' or 'chain_reaction'

class BioKernel:
    def __init__(self):
        self.miner = LiteratureMiningAgent()
        self.chemist = MoleculeEvolutionAgent()
        self.deputy = SafetyOfficerAgent()

    async def execute_chain_reaction(self, initial_query: str):
        """
        The 'Posse' Workflow:
        1. Miner finds target.
        2. Chemist designs drug.
        3. Deputy checks safety.
        """
        print(f"\n🤠 [Kernel] FORMING POSSE for: '{initial_query}'")
        results = []

        # --- Step 1: Literature Mining ---
        if bus: await bus.publish("step_start", {"step": "Mining"}, "Kernel")
        mining_result = self.miner.mine_for_targets(initial_query)
        town_square.write("latest_discovery", mining_result, "Miner")
        results.append(f"📰 Miner found: {mining_result.get('target_protein', 'Nothing')}")
        await asyncio.sleep(1) # Dramatic pause

        if mining_result.get("status") == "found":
            target = mining_result["target_protein"]
            
            # --- Step 2: Evolution (Triggered by Blackboard) ---
            if bus: await bus.publish("step_start", {"step": "Evolution"}, "Kernel")
            print(f"  ↪️ Triggering Chemist for target: {target}")
            drug_result = self.chemist.evolve(target)
            candidate = drug_result["top_candidate"]
            town_square.write("candidate_drug", candidate, "Chemist")
            results.append(f"🧬 Chemist designed: {candidate} (Score: {drug_result['score']:.2f})")
            await asyncio.sleep(1)

            # --- Step 3: Safety Check (Triggered by Blackboard) ---
            if bus: await bus.publish("step_start", {"step": "Safety"}, "Kernel")
            print(f"  ↪️ Triggering Deputy to inspect: {candidate}")
            safety_result = self.deputy.inspect_output(f"Proposed drug {candidate} for target {target}")
            town_square.write("safety_report", safety_result, "Deputy")
            
            status_icon = "✅" if safety_result["status"] == "approved" else "❌"
            results.append(f"{status_icon} Deputy says: {safety_result['status']}")

        return "\n".join(results)

kernel = BioKernel()

@app.post("/v1/frontier/chain_reaction")
async def run_chain(request: AgentRequest):
    report = await kernel.execute_chain_reaction(request.query)
    return {"report": report}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001) # Port 8001 for Frontier
