# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic.schema.messages import __all__

EXPECTED_ALL = [
    "AIMessage",
    "AIMessageChunk",
    "BaseMessage",
    "BaseMessageChunk",
    "ChatMessage",
    "ChatMessageChunk",
    "FunctionMessage",
    "FunctionMessageChunk",
    "HumanMessage",
    "HumanMessageChunk",
    "SystemMessage",
    "SystemMessageChunk",
    "ToolMessage",
    "ToolMessageChunk",
    "_message_from_dict",
    "_message_to_dict",
    "message_to_dict",
    "get_buffer_string",
    "merge_content",
    "messages_from_dict",
    "messages_to_dict",
    "AnyMessage",
]


def test_all_imports() -> None:
    assert set(__all__) == set(EXPECTED_ALL)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
