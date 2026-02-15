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
import json
import asyncio
import logging
from typing import Any, Dict, List
import os

# Add project root to path to import server
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))

from platform.core_kernel.server import kernel, AgentRequest
from platform.optimizer.usdl_transpiler import USDLTranspiler, USDLSpec, Provider

# Setup logging (stderr so it doesn't break JSON-RPC on stdout)
logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='[MCP] %(message)s')
logger = logging.getLogger("mcp_server")

class MCPServer:
    def __init__(self):
        self.kernel = kernel
        self.transpiler = USDLTranspiler()
        self.name = "corekernel-enterprise"
        self.version = "2026.3.0"

    async def handle_message(self, message: Dict[str, Any]):
        msg_type = message.get("method")
        msg_id = message.get("id")

        if msg_type == "initialize":
            return {
                "jsonrpc": "2.0",
                "id": msg_id,
                "result": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {
                        "tools": {},
                        "prompts": {}
                    },
                    "serverInfo": {
                        "name": self.name,
                        "version": self.version
                    }
                }
            }

        elif msg_type == "notifications/initialized":
            logger.info("Client initialized.")
            return None

        elif msg_type == "tools/list":
            tools = []
            # Tool 1: Execute Agent
            tools.append({
                "name": "run_agent",
                "description": "Executes an AI agent skill (e.g., Clinical Trial Matcher, CRISPR Designer).",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "query": {
                            "type": "string",
                            "description": "The natural language request for the agent."
                        },
                        "skill_id": {
                            "type": "string",
                            "description": "Optional: Specific skill ID if known (e.g., 'clinical-trial-matcher')."
                        }
                    },
                    "required": ["query"]
                }
            })
            # Tool 2: Transpile Skill (Universal Adapter)
            tools.append({
                "name": "transpile_skill",
                "description": "Converts a SKILL.md definition into an optimized prompt for OpenAI, Anthropic, or Gemini.",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "skill_name": {
                            "type": "string",
                            "description": "Name of the skill to transpile (e.g. 'clinical-trial-matcher')."
                        },
                        "target_provider": {
                            "type": "string",
                            "enum": ["openai", "anthropic", "gemini"],
                            "description": "The target LLM provider format."
                        }
                    },
                    "required": ["skill_name", "target_provider"]
                }
            })
            return {"jsonrpc": "2.0", "id": msg_id, "result": {"tools": tools}}

        elif msg_type == "tools/call":
            params = message.get("params", {})
            name = params.get("name")
            args = params.get("arguments", {})

            if name == "run_agent":
                query = args.get("query")
                skill_id = args.get("skill_id")
                
                # If skill_id not provided, let kernel route it
                if not skill_id:
                    req = AgentRequest(query=query)
                    skill_id = await self.kernel.route_request(req)

                result = await self.kernel.execute(skill_id, query, {})
                
                return {
                    "jsonrpc": "2.0", 
                    "id": msg_id, 
                    "result": {
                        "content": [{
                            "type": "text",
                            "text": result.response
                        }]
                    }
                }

            elif name == "transpile_skill":
                skill_name = args.get("skill_name")
                provider_str = args.get("target_provider")
                
                # Find skill
                skill_meta = self.kernel.skills_registry.get(skill_name)
                if not skill_meta:
                     return {"jsonrpc": "2.0", "id": msg_id, "error": {"code": -32001, "message": f"Skill '{skill_name}' not found."}}
                
                try:
                    # Load and Transpile
                    if skill_meta["type"] == "antigravity":
                        spec = USDLSpec.from_skill_md(skill_meta["content"])
                    else:
                        # Fallback for python agents - just use description
                        spec = USDLSpec(
                            name=skill_name, 
                            description="Legacy Python Agent", 
                            inputs=[], outputs=[], safety_checks=[], audit_policy=""
                        )

                    artifact = self.transpiler.compile(spec, Provider(provider_str))
                    
                    return {
                        "jsonrpc": "2.0",
                        "id": msg_id,
                        "result": {
                            "content": [{
                                "type": "text",
                                "text": json.dumps(artifact, indent=2)
                            }]
                        }
                    }
                except Exception as e:
                    return {"jsonrpc": "2.0", "id": msg_id, "error": {"code": -32002, "message": f"Transpilation failed: {str(e)}"}}

            else:
                return {"jsonrpc": "2.0", "id": msg_id, "error": {"code": -32601, "message": "Method not found"}}

        elif msg_type == "prompts/list":
            prompts = []
            # List all Antigravity skills as prompts
            for skill_id, meta in self.kernel.skills_registry.items():
                if meta.get("type") == "antigravity":
                    prompts.append({
                        "name": skill_id,
                        "description": meta.get("description"),
                        "arguments": []
                    })
            return {"jsonrpc": "2.0", "id": msg_id, "result": {"prompts": prompts}}

        elif msg_type == "prompts/get":
            params = message.get("params", {})
            name = params.get("name")
            skill_meta = self.kernel.skills_registry.get(name)

            if skill_meta and skill_meta.get("type") == "antigravity":
                return {
                    "jsonrpc": "2.0",
                    "id": msg_id,
                    "result": {
                        "messages": [{
                            "role": "user",
                            "content": {
                                "type": "text",
                                "text": skill_meta["content"]
                            }
                        }]
                    }
                }
            else:
                 return {"jsonrpc": "2.0", "id": msg_id, "error": {"code": -32602, "message": "Prompt not found"}}

        return None

    async def run(self):
        reader = asyncio.StreamReader()
        protocol = asyncio.StreamReaderProtocol(reader)
        await asyncio.get_running_loop().connect_read_pipe(lambda: protocol, sys.stdin)
        w_transport, w_protocol = await asyncio.get_running_loop().connect_write_pipe(asyncio.BaseProtocol, sys.stdout)

        logger.info("CoreKernel MCP Server running on stdio...")
        
        while True:
            try:
                line = await reader.readline()
                if not line:
                    break
                
                msg = json.loads(line)
                response = await self.handle_message(msg)
                
                if response:
                    sys.stdout.write(json.dumps(response) + "\n")
                    sys.stdout.flush()
            except Exception as e:
                logger.error(f"Error: {e}")

if __name__ == "__main__":
    server = MCPServer()
    try:
        asyncio.run(server.run())
    except KeyboardInterrupt:
        pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
