# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any

import pytest

from langchain_core.messages import AIMessage
from langchain_core.outputs import ChatGeneration


@pytest.mark.parametrize(
    "content",
    [
        "foo",
        ["foo"],
        [{"text": "foo", "type": "text"}],
        [
            {"tool_use": {}, "type": "tool_use"},
            {"text": "foo", "type": "text"},
            "bar",
        ],
    ],
)
def test_msg_with_text(content: str | list[str | dict[str, Any]]) -> None:
    expected = "foo"
    actual = ChatGeneration(message=AIMessage(content=content)).text
    assert actual == expected


@pytest.mark.parametrize("content", [[], [{"tool_use": {}, "type": "tool_use"}]])
def test_msg_no_text(content: str | list[str | dict[str, Any]]) -> None:
    expected = ""
    actual = ChatGeneration(message=AIMessage(content=content)).text
    assert actual == expected

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
