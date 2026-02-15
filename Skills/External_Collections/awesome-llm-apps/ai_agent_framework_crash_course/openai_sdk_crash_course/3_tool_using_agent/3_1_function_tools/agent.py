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
from .tools import add_numbers, multiply_numbers, get_weather, convert_temperature

# Create an agent with custom function tools
root_agent = Agent(
    name="Function Tools Agent",
    instructions="""
    You are a helpful assistant with access to various tools.
    
    Available tools:
    - add_numbers: Add two numbers together
    - multiply_numbers: Multiply two numbers together  
    - get_weather: Get weather information for a city
    - convert_temperature: Convert between Celsius and Fahrenheit
    
    When users ask for calculations or information:
    1. Use the appropriate tool for the task
    2. Explain what you're doing
    3. Show the result clearly
    
    Always use the provided tools rather than doing calculations yourself.
    """,
    tools=[add_numbers, multiply_numbers, get_weather, convert_temperature]
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
