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

app = FastAPI(title="CoreKernel Enterprise", version="2026.3.0-PRO")

class AgentRequest(BaseModel):
    query: str
    context: Optional[Dict[str, Any]] = {}
    model_preference: str = "auto" # 'auto', 'gemini', 'local', 'gpt4'

class AgentResponse(BaseModel):
    response: str
    tools_used: List[str]
    execution_time: float
    model_used: str

class CoreKernel:
    def __init__(self):
        self.skills_registry = {}
        self.medprompt = MedPromptEngine() if MedPromptEngine else None
        
        # Initialize WAL Providers
        self.providers = {}
        self._init_providers()

        self._discover_skills()
        self._discover_antigravity_skills()

    def _init_providers(self):
        # 1. Gemini 2.0 (Primary for speed/multimodal)
        gemini_key = os.getenv("GOOGLE_API_KEY")
        if gemini_key:
            self.providers["gemini"] = LLMFactory.create_provider("gemini", {"api_key": gemini_key, "model": "gemini-2.0-flash"})
            print("🔌 [CoreKernel] Loaded Provider: Gemini 2.0 Flash")

        # 2. Claude 3.7 (Primary for reasoning/thinking)
        anthropic_key = os.getenv("ANTHROPIC_API_KEY")
        if anthropic_key:
            self.providers["anthropic"] = LLMFactory.create_provider("anthropic", {"api_key": anthropic_key})
            print("🔌 [CoreKernel] Loaded Provider: Claude 3.7 Sonnet")

        # 3. Fallback
        self.providers["local"] = LLMFactory.create_provider("local", {})
        
        # Default LLM
        self.llm = self.providers.get("gemini") or self.providers.get("anthropic") or self.providers["local"]

    def _discover_skills(self):
        """
        Enterprise Mode: Dynamically scans the Skills/ directory for any file ending in '_agent.py'.
        """
        base_dir = "Skills"
        if not os.path.exists(base_dir):
            return

        print(f"🚀 [CoreKernel] Scanning {base_dir} for active agents...")
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
        print(f"🚀 [CoreKernel] Scanning for Antigravity SKILL.md files...")
        
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
            await bus.publish("kernel_routing", {"query": request.query}, "CoreKernel")

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
        
        # Determine Provider based on preference or complexity
        pref = context.get("model_preference", "auto")
        llm = self.providers.get(pref, self.llm)
        
        # Use Claude 3.7 for explicit 'thinking' requests
        if pref == "auto" and ("complex" in query.lower() or "design" in query.lower()):
            llm = self.providers.get("anthropic", self.llm)

        # 1. MedPrompt Injection for Clinical Queries
        if self.medprompt and ("patient" in query.lower() or "diagnosis" in query.lower()):
            print("  ⚕️ [CoreKernel] Injecting MedPrompt Strategy...")
            query = self.medprompt.generate_prompt(query)

        response_text = ""
        tools = []
        model_used = "unknown"

        # 2. Execution Logic
        skill_meta = self.skills_registry.get(skill_id)
        
        available_tools = [
            {"name": "literature_search", "description": "Search biomedical literature"},
            {"name": "mcp_github_tool", "description": "Interact with GitHub via MCP"},
            {"name": "mcp_filesystem_tool", "description": "Interact with local files via MCP"}
        ]

        print(f"  🧠 [CoreKernel] Reasoning via WAL (Provider: {getattr(llm, 'model', 'local')})...")
        
        # If it's an Antigravity skill, we pass its instructions to the Provider
        if skill_meta and skill_meta.get("type") == "antigravity":
            system_instr = f"Act as the following agent:\nName: {skill_meta['name']}\nDescription: {skill_meta['description']}\nInstructions: {skill_meta['content']}"
            
            # Pass reasoning budget if using Claude 3.7
            thinking_budget = context.get("thinking_budget", 4000) if "claude-3-7" in str(getattr(llm, 'model', '')) else 0

            req = LLMRequest(
                query=query, 
                system_instruction=system_instr,
                tools=available_tools,
                max_tokens=thinking_budget + 4000 if thinking_budget > 0 else 4000
            )
            llm_resp = await llm.generate(req)
            response_text = llm_resp.text
            model_used = llm_resp.model
            tools = ["antigravity_engine"]
        else:
            # Use general reasoning loop via WAL
            response_text = await llm.run_reasoning_loop(query, available_tools)
            model_used = str(getattr(llm, 'model', 'reasoning_engine'))
            tools = ["reasoning_engine"]

        # 3. Event Bus Notification
        if bus:
            await bus.publish("agent_completion", {
                "skill": skill_id, 
                "status": "success",
                "duration": time.time() - start_time
            }, "CoreKernel")

        return AgentResponse(
            response=response_text,
            tools_used=tools,
            execution_time=time.time() - start_time,
            model_used=model_used
        )

kernel = CoreKernel()

@app.post("/v1/agent/run", response_model=AgentResponse)
async def run_agent(request: AgentRequest):
    skill_id = await kernel.route_request(request)
    result = await kernel.execute(skill_id, request.query, request.context)
    return result

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
