# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest
from langchain_core.agents import AgentAction, AgentFinish
from langchain_core.exceptions import OutputParserException

from langchain_classic.agents.output_parsers.react_single_input import (
    ReActSingleInputOutputParser,
)


def test_action() -> None:
    """Test standard parsing of action/action input."""
    parser = ReActSingleInputOutputParser()
    _input = """Thought: agent thought here
Action: search
Action Input: what is the temperature in SF?"""
    output = parser.invoke(_input)
    expected_output = AgentAction(
        tool="search",
        tool_input="what is the temperature in SF?",
        log=_input,
    )
    assert output == expected_output


def test_finish() -> None:
    """Test standard parsing of agent finish."""
    parser = ReActSingleInputOutputParser()
    _input = """Thought: agent thought here
Final Answer: The temperature is 100"""
    output = parser.invoke(_input)
    expected_output = AgentFinish(
        return_values={"output": "The temperature is 100"},
        log=_input,
    )
    assert output == expected_output


def test_action_with_finish() -> None:
    """Test that if final thought is in action/action input, error is raised."""
    parser = ReActSingleInputOutputParser()
    _input = """Thought: agent thought here
Action: search Final Answer:
Action Input: what is the temperature in SF?"""
    with pytest.raises(OutputParserException):
        parser.invoke(_input)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
