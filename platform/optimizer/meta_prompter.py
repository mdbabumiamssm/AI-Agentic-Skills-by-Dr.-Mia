# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import re
from enum import Enum
from typing import Dict, List

# The Optimizer: Transforming Generic Prompts into SOTA Model-Specific Artifacts
# "One Prompt Source, Many Optimized Outputs"

class ModelTarget(Enum):
    CLAUDE = "claude"
    OPENAI = "openai"
    GEMINI = "gemini"

class PromptOptimizer:
    def __init__(self):
        pass

    def optimize(self, generic_prompt: str, target: ModelTarget) -> str:
        """
        Main entry point for optimization logic.
        """
        if target == ModelTarget.CLAUDE:
            return self._optimize_for_claude(generic_prompt)
        elif target == ModelTarget.OPENAI:
            return self._optimize_for_openai(generic_prompt)
        elif target == ModelTarget.GEMINI:
            return self._optimize_for_gemini(generic_prompt)
        return generic_prompt

    def _optimize_for_claude(self, prompt: str) -> str:
        """
        Applies Anthropic's specific best practices:
        - XML tagging for clear structure
        - Pre-fill injection (simulated here)
        - Thinking block encouragement
        """
        optimized = "<system_instructions>\n"
        optimized += "You are a highly capable biomedical AI assistant.\n"
        optimized += "Please think step-by-step inside <thinking> tags before answering.\n"
        optimized += "</system_instructions>\n\n"
        
        optimized += "<user_input>\n"
        optimized += prompt.strip()
        optimized += "\n</user_input>\n"
        
        return optimized

    def _optimize_for_openai(self, prompt: str) -> str:
        """
        Applies OpenAI's specific best practices:
        - Concise system messages
        - Explicit Markdown formatting instructions
        """
        # OpenAI prefers instruction at the top
        optimized = "### SYSTEM\n"
        optimized += "You are a helpful assistant for healthcare professionals.\n"
        optimized += "Provide concise, evidence-based answers.\n\n"
        
        optimized += "### USER\n"
        optimized += prompt.strip()
        
        return optimized

    def _optimize_for_gemini(self, prompt: str) -> str:
        """
        Applies Gemini's specific best practices:
        - Long context utilization
        - Structured input
        """
        # Gemini handles loose structure well, but benefits from explicit "Role"
        optimized = "Role: Expert Biomedical Researcher\n"
        optimized += "Task: Analyze the following input carefully.\n"
        optimized += "-" * 20 + "\n"
        optimized += prompt.strip() + "\n"
        optimized += "-" * 20 + "\n"
        optimized += "Requirement: Provide a detailed, multi-faceted response."
        
        return optimized

# Example Usage for Testing
if __name__ == "__main__":
    optimizer = PromptOptimizer()
    generic = "Analyze the side effects of Lisinopril for a 65-year-old male with diabetes."
    
    print("---> Claude Optimized ---")
    print(optimizer.optimize(generic, ModelTarget.CLAUDE))
    print("\n---> OpenAI Optimized ---")
    print(optimizer.optimize(generic, ModelTarget.OPENAI))
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
