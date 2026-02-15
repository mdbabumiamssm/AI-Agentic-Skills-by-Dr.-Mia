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

from langchain_classic.agents.mrkl.output_parser import (
    MISSING_ACTION_AFTER_THOUGHT_ERROR_MESSAGE,
    MISSING_ACTION_INPUT_AFTER_ACTION_ERROR_MESSAGE,
    MRKLOutputParser,
)

mrkl_output_parser = MRKLOutputParser()


def test_valid_action_and_action_input_parse() -> None:
    llm_output = """I can use the `foo` tool to achieve the goal.
    Action: foo
    Action Input: bar"""

    agent_action: AgentAction = mrkl_output_parser.parse(llm_output)  # type: ignore[assignment]
    assert agent_action.tool == "foo"
    assert agent_action.tool_input == "bar"


def test_valid_final_answer_parse() -> None:
    llm_output = """Final Answer: The best pizza to eat is margaritta """

    agent_finish: AgentFinish = mrkl_output_parser.parse(llm_output)  # type: ignore[assignment]
    assert (
        agent_finish.return_values.get("output")
        == "The best pizza to eat is margaritta"
    )


def test_missing_action() -> None:
    llm_output = """I can use the `foo` tool to achieve the goal."""

    with pytest.raises(OutputParserException) as exception_info:
        mrkl_output_parser.parse(llm_output)
    assert (
        exception_info.value.observation == MISSING_ACTION_AFTER_THOUGHT_ERROR_MESSAGE
    )


def test_missing_action_input() -> None:
    llm_output = """I can use the `foo` tool to achieve the goal.
    Action: foo"""

    with pytest.raises(OutputParserException) as exception_info:
        mrkl_output_parser.parse(llm_output)
    assert (
        exception_info.value.observation
        == MISSING_ACTION_INPUT_AFTER_ACTION_ERROR_MESSAGE
    )


def test_final_answer_before_parsable_action() -> None:
    llm_output = """Final Answer: The best pizza to eat is margaritta

        Action: foo
        Action Input: bar
        """
    agent_finish: AgentFinish = mrkl_output_parser.parse(llm_output)  # type: ignore[assignment]
    assert (
        agent_finish.return_values.get("output")
        == "The best pizza to eat is margaritta"
    )


def test_final_answer_after_parsable_action() -> None:
    llm_output = """
        Observation: I can use the `foo` tool to achieve the goal.
        Action: foo
        Action Input: bar
        Final Answer: The best pizza to eat is margaritta
        """
    with pytest.raises(OutputParserException) as exception_info:
        mrkl_output_parser.parse(llm_output)
    assert (
        "Parsing LLM output produced both a final answer and a parse-able action"
        in exception_info.value.args[0]
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
