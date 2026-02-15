# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Literal

from langchain_core.agents import AgentAction


def _escape(xml: str) -> str:
    """Replace XML tags with custom safe delimiters."""
    replacements = {
        "<tool>": "[[tool]]",
        "</tool>": "[[/tool]]",
        "<tool_input>": "[[tool_input]]",
        "</tool_input>": "[[/tool_input]]",
        "<observation>": "[[observation]]",
        "</observation>": "[[/observation]]",
    }
    for orig, repl in replacements.items():
        xml = xml.replace(orig, repl)
    return xml


def format_xml(
    intermediate_steps: list[tuple[AgentAction, str]],
    *,
    escape_format: Literal["minimal"] | None = "minimal",
) -> str:
    """Format the intermediate steps as XML.

    Args:
        intermediate_steps: The intermediate steps.
        escape_format: The escaping format to use. Currently only 'minimal' is
            supported, which replaces XML tags with custom delimiters to prevent
            conflicts.

    Returns:
        The intermediate steps as XML.
    """
    log = ""
    for action, observation in intermediate_steps:
        if escape_format == "minimal":
            # Escape XML tags in tool names and inputs using custom delimiters
            tool = _escape(action.tool)
            tool_input = _escape(str(action.tool_input))
            observation_ = _escape(str(observation))
        else:
            tool = action.tool
            tool_input = str(action.tool_input)
            observation_ = str(observation)
        log += (
            f"<tool>{tool}</tool><tool_input>{tool_input}"
            f"</tool_input><observation>{observation_}</observation>"
        )
    return log

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
