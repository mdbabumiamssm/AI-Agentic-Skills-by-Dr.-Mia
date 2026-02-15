# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

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

# Workflow Abstraction Layer (WAL) Imports
from platform.adapters.factory import LLMFactory
from platform.schema.io_types import LLMRequest, LLMResponse

app = FastAPI(title="BioKernel Enterprise", version="2026.3.0-PRO")

class AgentRequest(BaseModel):
    query: str
    context: Optional[Dict[str, Any]] = {}
    model_preference: str = "auto" # 'auto', 'gemini', 'local', 'gpt4'

class AgentResponse(BaseModel):
    response: str
    tools_used: List[str]
    execution_time: float
    model_used: str

class BioKernel:
    def __init__(self):
        self.skills_registry = {}
        self.medprompt = MedPromptEngine() if MedPromptEngine else None
        
        # Initialize WAL Provider (Default to Gemini, fallback to Local)
        # In production, this config would come from a YAML file or ENV
        self.provider_config = {
            "api_key": os.getenv("GOOGLE_API_KEY"),
            "model": "gemini-2.0-flash"
        }
        
        # Try to load primary provider, fallback to local if key missing
        if self.provider_config["api_key"]:
            print("🔌 [BioKernel] Loading Primary Provider: Gemini")
            self.llm = LLMFactory.create_provider("gemini", self.provider_config)
        else:
            print("🔌 [BioKernel] Loading Fallback Provider: Local (No API Key found)")
            self.llm = LLMFactory.create_provider("local", {})

        self._discover_skills()
        self._discover_antigravity_skills()

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
                    
                    self.skills_registry[skill_id] = {
                        "path": os.path.join(root, file),
                        "type": "python",
                        "name": skill_id
                    }
                    print(f"  + Registered (Python): {skill_id}")

    def _discover_antigravity_skills(self):
        """
        Antigravity Mode: Scans for SKILL.md files.
        """
        search_dirs = ["skill collections/Antigravity_Skills", "Skills"]
        print(f"🚀 [BioKernel] Scanning for Antigravity SKILL.md files...")
        
        for base_dir in search_dirs:
            if not os.path.exists(base_dir):
                continue
                
            for root, _, files in os.walk(base_dir):
                if "SKILL.md" in files:
                    file_path = os.path.join(root, "SKILL.md")
                    try:
                        with open(file_path, 'r') as f:
                            content = f.read()
                            # Simple frontmatter parser
                            if content.startswith("---"):
                                parts = content.split("---", 2)
                                if len(parts) >= 3:
                                    frontmatter_raw = parts[1]
                                    name = None
                                    description = "No description"
                                    
                                    for line in frontmatter_raw.splitlines():
                                        if line.strip().startswith("name:"):
                                            name = line.split(":", 1)[1].strip()
                                        elif line.strip().startswith("description:"):
                                            description = line.split(":", 1)[1].strip()
                                    
                                    if name:
                                        self.skills_registry[name] = {
                                            "path": file_path,
                                            "type": "antigravity",
                                            "description": description,
                                            "content": parts[2] # Markdown body
                                        }
                                        print(f"  + Registered (Antigravity): {name}")
                    except Exception as e:
                        print(f"  ! Failed to load {file_path}: {e}")

    async def route_request(self, request: AgentRequest) -> str:
        """
        Intelligent Router with Event Bus logging.
        """
        query_lower = request.query.lower()
        
        if bus:
            await bus.publish("kernel_routing", {"query": request.query}, "BioKernel")

        # 1. Check Antigravity Skills first (Direct Match)
        for skill_id, meta in self.skills_registry.items():
            if meta["type"] == "antigravity":
                # Simple keyword matching based on description or name
                keywords = skill_id.split("-") + meta["description"].lower().split()
                if any(k in query_lower for k in keywords if len(k) > 3):
                    return skill_id

        # 2. Fallback to Python Agents
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

        response_text = ""
        tools = []
        model_used = "unknown"

        # 2. Execution Logic
        skill_meta = self.skills_registry.get(skill_id)
        
        # Define available tools (in a real system, these would be dynamic)
        available_tools = [
            {"name": "literature_search", "description": "Search biomedical literature"},
            {"name": "clinical_trial_match", "description": "Find matching clinical trials"},
            {"name": "molecule_generator", "description": "Generate small molecules for a target"}
        ]

        print(f"  🧠 [BioKernel] Reasoning via WAL...")
        
        # If it's an Antigravity skill, we pass its instructions to the Provider
        if skill_meta and skill_meta.get("type") == "antigravity":
            system_instr = f"Act as the following agent:\nName: {skill_meta['name']}\nDescription: {skill_meta['description']}\nInstructions: {skill_meta['content']}"
            
            req = LLMRequest(
                query=query, 
                system_instruction=system_instr,
                tools=available_tools # WAL supports tools
            )
            llm_resp = await self.llm.generate(req)
            response_text = llm_resp.text
            model_used = llm_resp.model
            tools = ["antigravity_engine"]
        else:
            # Use general reasoning loop via WAL
            # The provider (Gemini/Local) handles the specific prompt engineering
            response_text = await self.llm.run_reasoning_loop(query, available_tools)
            model_used = "reasoning_engine"
            tools = ["reasoning_engine"]

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
            model_used=model_used
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

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
