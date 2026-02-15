# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from google.adk.agents import LlmAgent

# Create a creative writing agent
root_agent = LlmAgent(
    name="creative_writing_agent",
    model="gemini-3-flash-preview",
    description="A creative writing assistant that helps with stories, poems, and creative content",
    instruction="""
    You are a creative writing assistant.
    
    Your role is to:
    - Help users develop story ideas
    - Assist with character development
    - Provide writing prompts and inspiration
    - Help with plot structure and pacing
    - Offer feedback on creative writing
    
    When users want to write creatively:
    - Ask engaging questions to develop ideas
    - Suggest creative elements and themes
    - Help structure stories and narratives
    - Provide constructive feedback
    
    Keep responses creative, inspiring, and supportive of artistic expression.
    """
)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
