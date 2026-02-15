# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.agents import AgentAction
from langchain_core.messages import AIMessage, BaseMessage, HumanMessage


def format_log_to_messages(
    intermediate_steps: list[tuple[AgentAction, str]],
    template_tool_response: str = "{observation}",
) -> list[BaseMessage]:
    """Construct the scratchpad that lets the agent continue its thought process.

    Args:
        intermediate_steps: List of tuples of AgentAction and observation strings.
        template_tool_response: Template to format the observation with.
            Defaults to `"{observation}"`.

    Returns:
        The scratchpad.
    """
    thoughts: list[BaseMessage] = []
    for action, observation in intermediate_steps:
        thoughts.append(AIMessage(content=action.log))
        human_message = HumanMessage(
            content=template_tool_response.format(observation=observation),
        )
        thoughts.append(human_message)
    return thoughts

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
