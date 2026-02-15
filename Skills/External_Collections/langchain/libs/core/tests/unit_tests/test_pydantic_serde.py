# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test pydantic SerDe.

A set of tests that verifies that Union discrimination works correctly with
the various pydantic base models.

These tests can uncover issues that will also arise during regular instantiation
of the models (i.e., not necessarily from loading or dumping JSON).
"""

import pytest
from pydantic import RootModel, ValidationError

from langchain_core.messages import (
    AIMessage,
    AIMessageChunk,
    AnyMessage,
    ChatMessage,
    ChatMessageChunk,
    FunctionMessage,
    FunctionMessageChunk,
    HumanMessage,
    HumanMessageChunk,
    SystemMessage,
    SystemMessageChunk,
)


def test_serde_any_message() -> None:
    """Test AnyMessage() serder."""
    lc_objects = [
        HumanMessage(content="human"),
        HumanMessageChunk(content="human"),
        AIMessage(content="ai"),
        AIMessageChunk(content="ai"),
        SystemMessage(content="sys"),
        SystemMessageChunk(content="sys"),
        FunctionMessage(
            name="func",
            content="func",
        ),
        FunctionMessageChunk(
            name="func",
            content="func",
        ),
        ChatMessage(
            role="human",
            content="human",
        ),
        ChatMessageChunk(
            role="human",
            content="human",
        ),
    ]

    model = RootModel[AnyMessage]

    for lc_object in lc_objects:
        d = lc_object.model_dump()
        assert "type" in d, f"Missing key `type` for {type(lc_object)}"
        obj1 = model.model_validate(d)
        assert type(obj1.root) is type(lc_object), f"failed for {type(lc_object)}"

    with pytest.raises((TypeError, ValidationError)):
        # Make sure that specifically validation error is raised
        model.model_validate({})

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
