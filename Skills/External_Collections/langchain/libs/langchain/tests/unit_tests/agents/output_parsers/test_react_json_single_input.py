# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.agents import AgentAction, AgentFinish

from langchain_classic.agents.output_parsers.react_json_single_input import (
    ReActJsonSingleInputOutputParser,
)


def test_action() -> None:
    """Test standard parsing of action/action input."""
    parser = ReActJsonSingleInputOutputParser()
    _input = """Thought: agent thought here
```
{
    "action": "search",
    "action_input": "what is the temperature in SF?"
}
```
"""
    output = parser.invoke(_input)
    expected_output = AgentAction(
        tool="search",
        tool_input="what is the temperature in SF?",
        log=_input,
    )
    assert output == expected_output


def test_finish() -> None:
    """Test standard parsing of agent finish."""
    parser = ReActJsonSingleInputOutputParser()
    _input = """Thought: agent thought here
Final Answer: The temperature is 100"""
    output = parser.invoke(_input)
    expected_output = AgentFinish(
        return_values={"output": "The temperature is 100"},
        log=_input,
    )
    assert output == expected_output

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
