# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from agents import Agent
from agents.tools import WebSearchTool, CodeInterpreterTool

# Create an agent with built-in OpenAI tools
root_agent = Agent(
    name="Built-in Tools Agent",
    instructions="""
    You are a research and computation assistant with access to powerful built-in tools.
    
    Available tools:
    - WebSearchTool: Search the web for current information
    - CodeInterpreterTool: Execute Python code safely
    
    You can help with:
    - Finding current information and news
    - Performing complex calculations
    - Data analysis and visualization
    - Mathematical computations
    
    When users request information or calculations:
    1. Use web search for current information
    2. Use code execution for computations and analysis
    3. Provide clear explanations of results
    """,
    tools=[WebSearchTool(), CodeInterpreterTool()]
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
