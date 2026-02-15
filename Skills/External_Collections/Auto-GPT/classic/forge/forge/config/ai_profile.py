# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from pydantic import BaseModel, Field

DEFAULT_AI_NAME = "AutoGPT"
DEFAULT_AI_ROLE = (
    "a seasoned digital assistant: "
    "capable, intelligent, considerate and assertive. "
    "You have extensive research and development skills, and you don't shy "
    "away from writing some code to solve a problem. "
    "You are pragmatic and make the most out of the tools available to you."
)


class AIProfile(BaseModel):
    """
    Object to hold the AI's personality.

    Attributes:
        ai_name (str): The name of the AI.
        ai_role (str): The description of the AI's role.
        ai_goals (list): The list of objectives the AI is supposed to complete.
        api_budget (float): The maximum dollar value for API calls (0.0 means infinite)
    """

    ai_name: str = DEFAULT_AI_NAME
    ai_role: str = DEFAULT_AI_ROLE
    """`ai_role` should fit in the following format: `You are {ai_name}, {ai_role}`"""
    ai_goals: list[str] = Field(default_factory=list[str])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
