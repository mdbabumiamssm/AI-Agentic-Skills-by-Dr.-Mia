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

from langchain_classic.base_memory import BaseMemory
from langchain_classic.chains.conversation.memory import (
    ConversationBufferMemory,
    ConversationBufferWindowMemory,
    ConversationSummaryMemory,
)
from langchain_classic.memory import ReadOnlySharedMemory, SimpleMemory
from tests.unit_tests.llms.fake_llm import FakeLLM


def test_simple_memory() -> None:
    """Test SimpleMemory."""
    memory = SimpleMemory(memories={"baz": "foo"})

    output = memory.load_memory_variables({})

    assert output == {"baz": "foo"}
    assert memory.memory_variables == ["baz"]


@pytest.mark.parametrize(
    "memory",
    [
        ConversationBufferMemory(memory_key="baz"),
        ConversationSummaryMemory(llm=FakeLLM(), memory_key="baz"),
        ConversationBufferWindowMemory(memory_key="baz"),
    ],
)
def test_readonly_memory(memory: BaseMemory) -> None:
    read_only_memory = ReadOnlySharedMemory(memory=memory)
    memory.save_context({"input": "bar"}, {"output": "foo"})

    assert read_only_memory.load_memory_variables({}) == memory.load_memory_variables(
        {},
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
