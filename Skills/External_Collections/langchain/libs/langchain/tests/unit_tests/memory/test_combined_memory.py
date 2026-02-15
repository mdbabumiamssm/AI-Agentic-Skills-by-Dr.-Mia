# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test for CombinedMemory class."""

import re

import pytest

from langchain_classic.memory import CombinedMemory, ConversationBufferMemory


@pytest.fixture
def example_memory() -> list[ConversationBufferMemory]:
    example_1 = ConversationBufferMemory(memory_key="foo")
    example_2 = ConversationBufferMemory(memory_key="bar")
    example_3 = ConversationBufferMemory(memory_key="bar")
    return [example_1, example_2, example_3]


def test_basic_functionality(example_memory: list[ConversationBufferMemory]) -> None:
    """Test basic functionality of methods exposed by class."""
    combined_memory = CombinedMemory(memories=[example_memory[0], example_memory[1]])
    assert combined_memory.memory_variables == ["foo", "bar"]
    assert combined_memory.load_memory_variables({}) == {"foo": "", "bar": ""}
    combined_memory.save_context(
        {"input": "Hello there"},
        {"output": "Hello, how can I help you?"},
    )
    assert combined_memory.load_memory_variables({}) == {
        "foo": "Human: Hello there\nAI: Hello, how can I help you?",
        "bar": "Human: Hello there\nAI: Hello, how can I help you?",
    }
    combined_memory.clear()
    assert combined_memory.load_memory_variables({}) == {"foo": "", "bar": ""}


def test_repeated_memory_var(example_memory: list[ConversationBufferMemory]) -> None:
    """Test raising error when repeated memory variables found."""
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Value error, The same variables {'bar'} are found in "
            "multiplememory object, which is not allowed by CombinedMemory."
        ),
    ):
        CombinedMemory(memories=[example_memory[1], example_memory[2]])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
