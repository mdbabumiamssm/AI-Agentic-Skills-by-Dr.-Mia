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

from langchain_classic.agents.output_parsers.self_ask import SelfAskOutputParser


def test_follow_up() -> None:
    """Test follow up parsing."""
    parser = SelfAskOutputParser()
    _input = "Follow up: what is two + 2"
    output = parser.invoke(_input)
    expected_output = AgentAction(
        tool="Intermediate Answer",
        tool_input="what is two + 2",
        log=_input,
    )
    assert output == expected_output
    # Test that also handles one word by default
    _input = "Followup: what is two + 2"
    output = parser.invoke(_input)
    expected_output = AgentAction(
        tool="Intermediate Answer",
        tool_input="what is two + 2",
        log=_input,
    )
    assert output == expected_output


def test_follow_up_custom() -> None:
    """Test follow up parsing for custom followups."""
    parser = SelfAskOutputParser(followups=("Now:",))
    _input = "Now: what is two + 2"
    output = parser.invoke(_input)
    expected_output = AgentAction(
        tool="Intermediate Answer",
        tool_input="what is two + 2",
        log=_input,
    )
    assert output == expected_output


def test_finish() -> None:
    """Test standard finish."""
    parser = SelfAskOutputParser()
    _input = "So the final answer is: 4"
    output = parser.invoke(_input)
    expected_output = AgentFinish(return_values={"output": "4"}, log=_input)
    assert output == expected_output


def test_finish_custom() -> None:
    """Test custom finish."""
    parser = SelfAskOutputParser(finish_string="Finally: ")
    _input = "Finally: 4"
    output = parser.invoke(_input)
    expected_output = AgentFinish(return_values={"output": "4"}, log=_input)
    assert output == expected_output

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
