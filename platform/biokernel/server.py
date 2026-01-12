import os
import json
import asyncio
import time
import importlib.util
from typing import Dict, Any, List, Optional
from fastapi import FastAPI
from pydantic import BaseModel

# Import Core Systems
try:
    from Skills.Computer_Science.Distributed_Systems.event_bus import bus
    from Skills.LLM_Research.Prompt_Engineering.medprompt import MedPromptEngine
except ImportError:
    # Fallback for when running in isolation/dev mode
    bus = None
    MedPromptEngine = None

app = FastAPI(title="BioKernel Enterprise", version="2026.3.0-PRO")

class AgentRequest(BaseModel):
    query: str
    context: Optional[Dict[str, Any]] = {}
    model_preference: str = "auto"

class AgentResponse(BaseModel):
    response: str
    tools_used: List[str]
    execution_time: float
    model_used: str

class BioKernel:
    def __init__(self):
        self.skills_registry = {}
        self.medprompt = MedPromptEngine() if MedPromptEngine else None
        self._discover_skills()

    def _discover_skills(self):
        """
        Enterprise Mode: Dynamically scans the Skills/ directory for any file ending in '_agent.py'.
        """
        base_dir = "Skills"
        if not os.path.exists(base_dir):
            return

        print(f"🚀 [BioKernel] Scanning {base_dir} for active agents...")
        for root, _, files in os.walk(base_dir):
            for file in files:
                if file.endswith("_agent.py") or file == "agent.py":
                    # Construct a unique ID like 'clinical_prior_auth'
                    rel_path = os.path.relpath(os.path.join(root, file), base_dir)
                    skill_id = rel_path.replace("/", "_").replace(".py", "").lower()
                    
                    self.skills_registry[skill_id] = os.path.join(root, file)
                    print(f"  + Registered: {skill_id}")

    async def route_request(self, request: AgentRequest) -> str:
        """
        Intelligent Router with Event Bus logging.
        """
        query_lower = request.query.lower()
        
        if bus:
            await bus.publish("kernel_routing", {"query": request.query}, "BioKernel")

        # Routing Logic
        if "insurance" in query_lower:
            return "clinical_prior_authorization_agent" # Matches file path structure
        elif "heart" in query_lower:
            return "consumer_health_wearable_analysis_agent"
        elif "molecule" in query_lower:
            return "drug_discovery_molecule_design_evolution_agent" # New agent
        elif "trial" in query_lower:
            return "clinical_clinical_trials_recruitment_agent" # New agent
        
        return "general_assistant"

    async def execute(self, skill_id: str, query: str, context: Dict) -> AgentResponse:
        start_time = time.time()
        
        # 1. MedPrompt Injection for Clinical Queries
        if self.medprompt and ("patient" in query.lower() or "diagnosis" in query.lower()):
            print("  ⚕️ [BioKernel] Injecting MedPrompt Strategy...")
            query = self.medprompt.generate_prompt(query)

        # 2. Execution (Mocked for now, but wired for expansion)
        # In a real system, we'd dynamically import and run the agent class here.
        # For this demo, we simulate the specific responses of the new agents.
        
        response_text = ""
        tools = []
        model = "claude-3-5-sonnet"

        if "evolution" in skill_id:
            response_text = "🧬 Evolution Complete: Generated 50 variants. Top candidate (Score: 0.98) binds to Target X."
            tools = ["genetic_algorithm", "docking_oracle"]
        elif "recruitment" in skill_id:
            response_text = "📋 Trial Match: Found 12 eligible patients from EHR. 3 High-priority matches."
            tools = ["vector_db_search", "fhir_parser"]
        elif "prior_auth" in skill_id:
            response_text = "✅ APPROVED: Meets criteria 2.1 (Conservative Therapy Failed)."
            tools = ["policy_engine"]
        else:
            response_text = f"Processed request using {skill_id}: {query[:50]}..."
            model = "gpt-4o"

        # 3. Event Bus Notification
        if bus:
            await bus.publish("agent_completion", {
                "skill": skill_id, 
                "status": "success",
                "duration": time.time() - start_time
            }, "BioKernel")

        return AgentResponse(
            response=response_text,
            tools_used=tools,
            execution_time=time.time() - start_time,
            model_used=model
        )

kernel = BioKernel()

@app.post("/v1/agent/run", response_model=AgentResponse)
async def run_agent(request: AgentRequest):
    skill_id = await kernel.route_request(request)
    result = await kernel.execute(skill_id, request.query, request.context)
    return result

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
