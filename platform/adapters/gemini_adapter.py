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
import time
from typing import List, Dict, Any, Optional
from platform.interface.llm_provider import LLMProvider
from platform.schema.io_types import LLMRequest, LLMResponse, ToolCall

try:
    import google.generativeai as genai
    HAS_GEMINI = True
except ImportError:
    HAS_GEMINI = False

class GeminiAdapter(LLMProvider):
    """
    Implements the LLMProvider interface for Google's Gemini models.
    """
    def __init__(self):
        self.model = None
        self.api_key = None
        self.available = False
        self.config = {}

    def initialize(self, config: Dict[str, Any]) -> bool:
        """
        Loads the API key and configures the model.
        """
        self.config = config
        self.api_key = config.get("api_key") or os.getenv("GOOGLE_API_KEY")
        
        if HAS_GEMINI and self.api_key:
            genai.configure(api_key=self.api_key)
            model_name = config.get("model", "gemini-2.0-flash")
            self.model = genai.GenerativeModel(model_name)
            self.available = True
            print(f"✅ [GeminiAdapter] Initialized with model: {model_name}")
            return True
        else:
            print("⚠️ [GeminiAdapter] Failed to initialize (missing key or dependency).")
            self.available = False
            return False

    async def generate(self, request: LLMRequest) -> LLMResponse:
        start_time = time.time()
        
        if not self.available:
            return LLMResponse(
                text="Error: Gemini adapter not available.",
                provider="gemini",
                model="unknown",
                finish_reason="error"
            )

        try:
            # Handle System Instruction
            model = self.model
            if request.system_instruction:
                model = genai.GenerativeModel(
                    self.config.get("model", "gemini-2.0-flash"),
                    system_instruction=request.system_instruction
                )

            # Generate
            response = await model.generate_content_async(
                request.query,
                generation_config=genai.types.GenerationConfig(
                    temperature=request.temperature,
                    max_output_tokens=request.max_tokens,
                    stop_sequences=request.stop_sequences
                )
            )
            
            latency = (time.time() - start_time) * 1000
            
            # Basic usage tracking (mocked as Gemini API doesn't always return this cleanly in all versions)
            usage = {"prompt_tokens": len(request.query) // 4, "completion_tokens": len(response.text) // 4}

            return LLMResponse(
                text=response.text,
                provider="gemini",
                model=self.config.get("model", "gemini-2.0-flash"),
                latency_ms=latency,
                usage=usage
            )

        except Exception as e:
            return LLMResponse(
                text=f"Error generating content: {str(e)}",
                provider="gemini",
                model="unknown",
                finish_reason="error"
            )

    async def run_reasoning_loop(self, goal: str, tools: List[Dict]) -> str:
        """
        Standard ReAct loop implementation using Gemini.
        """
        if not self.available:
            return "Simulation: Adapter unavailable."

        tool_desc = "\n".join([f"- {t['name']}: {t['description']}" for t in tools])
        prompt = f"""
You are an intelligent Biomedical Agent.
Goal: {goal}

Available Tools:
{tool_desc}

Think step-by-step. 
1. Analyze the goal.
2. Decide if you need a tool.
3. If yes, output METHOD: <tool_name> INPUT: <input>.
4. If you have the answer, output FINAL ANSWER: <answer>.
"""
        # Create a simplified request for the loop
        req = LLMRequest(query=prompt, temperature=0.2)
        resp = await self.generate(req)
        return resp.text
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
