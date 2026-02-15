# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from agents import Agent, Runner, RunConfig

# Create an agent for demonstrating run configuration
root_agent = Agent(
    name="Configuration Demo Agent",
    instructions="You are a helpful assistant that demonstrates run configuration options."
)

# Example 1: Basic run configuration with model settings
async def model_config_example():
    """Demonstrates run configuration with model overrides and settings"""
    
    run_config = RunConfig(
        model="gpt-4o",  # Override agent's default model
        model_settings={
            "temperature": 0.1,  # Low temperature for consistent responses
            "top_p": 0.9
        },
        max_turns=5,  # Limit conversation turns
        workflow_name="demo_workflow",  # For tracing
        trace_metadata={"experiment": "config_demo"}
    )
    
    result = await Runner.run(
        root_agent, 
        "Explain the weather in exactly 3 sentences.",
        run_config=run_config
    )
    
    return result.final_output

# Example 2: Run configuration with tracing settings
async def tracing_config_example():
    """Demonstrates run configuration with tracing options"""
    
    run_config = RunConfig(
        tracing_disabled=False,  # Enable tracing
        trace_include_sensitive_data=False,  # Exclude sensitive data
        workflow_name="production_workflow",
        group_id="user_session_456",  # Link multiple runs
        trace_metadata={
            "user_id": "user_123",
            "feature": "chat_assistance"
        }
    )
    
    result = await Runner.run(
        root_agent,
        "What are the benefits of structured logging?",
        run_config=run_config
    )
    
    return result.final_output

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
