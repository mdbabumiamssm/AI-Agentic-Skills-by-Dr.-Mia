# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test in memory indexer."""

from collections.abc import AsyncGenerator, Generator

import pytest
from langchain_tests.integration_tests.indexer import (
    AsyncDocumentIndexTestSuite,
    DocumentIndexerTestSuite,
)
from typing_extensions import override

from langchain_core.documents import Document
from langchain_core.indexing.base import DocumentIndex
from langchain_core.indexing.in_memory import (
    InMemoryDocumentIndex,
)


class TestDocumentIndexerTestSuite(DocumentIndexerTestSuite):
    @pytest.fixture
    @override
    def index(self) -> Generator[DocumentIndex, None, None]:
        yield InMemoryDocumentIndex()  # noqa: PT022


class TestAsyncDocumentIndexerTestSuite(AsyncDocumentIndexTestSuite):
    # Something funky is going on with mypy and async pytest fixture
    @pytest.fixture
    @override
    async def index(self) -> AsyncGenerator[DocumentIndex, None]:
        yield InMemoryDocumentIndex()  # noqa: PT022


def test_sync_retriever() -> None:
    index = InMemoryDocumentIndex()
    documents = [
        Document(id="1", page_content="hello world"),
        Document(id="2", page_content="goodbye cat"),
    ]
    index.upsert(documents)
    assert index.invoke("hello") == [documents[0], documents[1]]
    assert index.invoke("cat") == [documents[1], documents[0]]


async def test_async_retriever() -> None:
    index = InMemoryDocumentIndex()
    documents = [
        Document(id="1", page_content="hello world"),
        Document(id="2", page_content="goodbye cat"),
    ]
    await index.aupsert(documents)
    assert (await index.ainvoke("hello")) == [documents[0], documents[1]]
    assert (await index.ainvoke("cat")) == [documents[1], documents[0]]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
