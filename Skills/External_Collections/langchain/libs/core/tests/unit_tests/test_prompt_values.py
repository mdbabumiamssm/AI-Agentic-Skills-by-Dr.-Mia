# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.messages import (
    AIMessage,
    AIMessageChunk,
    HumanMessage,
    HumanMessageChunk,
    SystemMessage,
    SystemMessageChunk,
    ToolMessage,
    ToolMessageChunk,
)
from langchain_core.prompt_values import ChatPromptValueConcrete


def test_chat_prompt_value_concrete() -> None:
    messages: list = [
        AIMessage("foo"),
        HumanMessage("foo"),
        SystemMessage("foo"),
        ToolMessage("foo", tool_call_id="foo"),
        AIMessageChunk(content="foo"),
        HumanMessageChunk(content="foo"),
        SystemMessageChunk(content="foo"),
        ToolMessageChunk(content="foo", tool_call_id="foo"),
    ]
    assert ChatPromptValueConcrete(messages=messages).messages == messages

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
