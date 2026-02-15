# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test the initialize module."""

from langchain_core.tools import tool

from langchain_classic.agents.agent_types import AgentType
from langchain_classic.agents.initialize import initialize_agent
from tests.unit_tests.llms.fake_llm import FakeLLM


@tool
def my_tool(query: str) -> str:  # noqa: ARG001
    """A fake tool."""
    return "fake tool"


def test_initialize_agent_with_str_agent_type() -> None:
    """Test initialize_agent with a string."""
    fake_llm = FakeLLM()
    agent_executor = initialize_agent(
        [my_tool],
        fake_llm,
        "zero-shot-react-description",  # type: ignore[arg-type]
    )
    assert (
        agent_executor._action_agent._agent_type
        == AgentType.ZERO_SHOT_REACT_DESCRIPTION
    )
    assert isinstance(agent_executor.tags, list)
    assert "zero-shot-react-description" in agent_executor.tags

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
