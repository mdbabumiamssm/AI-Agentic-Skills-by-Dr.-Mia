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
from langchain_core.documents import Document

from langchain_classic.retrievers.multi_query import (
    LineListOutputParser,
    _unique_documents,
)


@pytest.mark.parametrize(
    ("documents", "expected"),
    [
        ([], []),
        ([Document(page_content="foo")], [Document(page_content="foo")]),
        ([Document(page_content="foo")] * 2, [Document(page_content="foo")]),
        (
            [Document(page_content="foo", metadata={"bar": "baz"})] * 2,
            [Document(page_content="foo", metadata={"bar": "baz"})],
        ),
        (
            [Document(page_content="foo", metadata={"bar": [1, 2]})] * 2,
            [Document(page_content="foo", metadata={"bar": [1, 2]})],
        ),
        (
            [Document(page_content="foo", metadata={"bar": {1, 2}})] * 2,
            [Document(page_content="foo", metadata={"bar": {1, 2}})],
        ),
        (
            [
                Document(page_content="foo", metadata={"bar": [1, 2]}),
                Document(page_content="foo", metadata={"bar": [2, 1]}),
            ],
            [
                Document(page_content="foo", metadata={"bar": [1, 2]}),
                Document(page_content="foo", metadata={"bar": [2, 1]}),
            ],
        ),
    ],
)
def test__unique_documents(documents: list[Document], expected: list[Document]) -> None:
    assert _unique_documents(documents) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("foo\nbar\nbaz", ["foo", "bar", "baz"]),
        ("foo\nbar\nbaz\n", ["foo", "bar", "baz"]),
        ("foo\n\nbar", ["foo", "bar"]),
    ],
)
def test_line_list_output_parser(text: str, expected: list[str]) -> None:
    parser = LineListOutputParser()
    assert parser.parse(text) == expected

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
