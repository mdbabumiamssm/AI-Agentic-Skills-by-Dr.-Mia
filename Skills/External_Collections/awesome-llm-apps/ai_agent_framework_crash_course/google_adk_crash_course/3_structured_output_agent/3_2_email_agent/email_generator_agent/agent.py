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
from pydantic import BaseModel, Field

class EmailContent(BaseModel):
    """Schema for email content with subject and body."""
    
    subject: str = Field(
        description="The subject line of the email. Should be concise and descriptive."
    )
    body: str = Field(
        description="The main content of the email. Should be well-formatted with proper greeting, paragraphs, and signature."
    )

root_agent = LlmAgent(
    name="email_generator_agent",
    model="gemini-3-flash-preview",
    description="Professional email generator that creates structured email content",
    instruction="""
    You are a professional email writing assistant. 
    
    IMPORTANT: Your response must be a JSON object with exactly these fields:
    - "subject": A concise, relevant subject line
    - "body": Well-formatted email content with greeting, main content, and closing
    
    Format your response as valid JSON only.
    """,
    output_schema=EmailContent,  # This is where the magic happens
    output_key="generated_email"  
)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
