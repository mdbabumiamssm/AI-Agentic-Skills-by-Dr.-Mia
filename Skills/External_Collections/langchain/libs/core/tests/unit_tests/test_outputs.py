# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.messages import HumanMessageChunk
from langchain_core.outputs import ChatGenerationChunk, GenerationChunk


def test_generation_chunk() -> None:
    assert GenerationChunk(text="Hello, ") + GenerationChunk(
        text="world!"
    ) == GenerationChunk(text="Hello, world!"), (
        "GenerationChunk + GenerationChunk should be a GenerationChunk"
    )

    assert GenerationChunk(text="Hello, ") + GenerationChunk(
        text="world!", generation_info={"foo": "bar"}
    ) == GenerationChunk(text="Hello, world!", generation_info={"foo": "bar"}), (
        "GenerationChunk + GenerationChunk should be a GenerationChunk "
        "with merged generation_info"
    )

    assert GenerationChunk(text="Hello, ") + GenerationChunk(
        text="world!", generation_info={"foo": "bar"}
    ) + GenerationChunk(text="!", generation_info={"baz": "foo"}) == GenerationChunk(
        text="Hello, world!!", generation_info={"foo": "bar", "baz": "foo"}
    ), (
        "GenerationChunk + GenerationChunk should be a GenerationChunk "
        "with merged generation_info"
    )


def test_chat_generation_chunk() -> None:
    assert ChatGenerationChunk(
        message=HumanMessageChunk(content="Hello, ")
    ) + ChatGenerationChunk(
        message=HumanMessageChunk(content="world!")
    ) == ChatGenerationChunk(message=HumanMessageChunk(content="Hello, world!")), (
        "ChatGenerationChunk + ChatGenerationChunk should be a ChatGenerationChunk"
    )

    assert ChatGenerationChunk(
        message=HumanMessageChunk(content="Hello, ")
    ) + ChatGenerationChunk(
        message=HumanMessageChunk(content="world!"), generation_info={"foo": "bar"}
    ) == ChatGenerationChunk(
        message=HumanMessageChunk(content="Hello, world!"),
        generation_info={"foo": "bar"},
    ), (
        "GenerationChunk + GenerationChunk should be a GenerationChunk "
        "with merged generation_info"
    )

    assert ChatGenerationChunk(
        message=HumanMessageChunk(content="Hello, ")
    ) + ChatGenerationChunk(
        message=HumanMessageChunk(content="world!"), generation_info={"foo": "bar"}
    ) + ChatGenerationChunk(
        message=HumanMessageChunk(content="!"), generation_info={"baz": "foo"}
    ) == ChatGenerationChunk(
        message=HumanMessageChunk(content="Hello, world!!"),
        generation_info={"foo": "bar", "baz": "foo"},
    ), (
        "GenerationChunk + GenerationChunk should be a GenerationChunk "
        "with merged generation_info"
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
