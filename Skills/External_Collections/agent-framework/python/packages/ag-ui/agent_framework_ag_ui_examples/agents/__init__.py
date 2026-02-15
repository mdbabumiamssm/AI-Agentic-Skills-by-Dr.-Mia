# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

"""Example agents for AG-UI demonstration."""

from .document_writer_agent import document_writer_agent
from .human_in_the_loop_agent import human_in_the_loop_agent
from .recipe_agent import recipe_agent
from .research_assistant_agent import research_assistant_agent
from .simple_agent import simple_agent
from .task_planner_agent import task_planner_agent
from .task_steps_agent import task_steps_agent_wrapped
from .ui_generator_agent import ui_generator_agent
from .weather_agent import weather_agent

__all__ = [
    "document_writer_agent",
    "human_in_the_loop_agent",
    "recipe_agent",
    "research_assistant_agent",
    "simple_agent",
    "task_planner_agent",
    "task_steps_agent_wrapped",
    "ui_generator_agent",
    "weather_agent",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
