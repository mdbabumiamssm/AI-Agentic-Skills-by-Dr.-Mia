# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import tempfile
from collections.abc import Generator
from typing import cast

import pytest
from langchain_core.documents import Document

from langchain_classic.storage._lc_store import create_kv_docstore, create_lc_store
from langchain_classic.storage.file_system import LocalFileStore


@pytest.fixture
def file_store() -> Generator[LocalFileStore, None, None]:
    # Create a temporary directory for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Instantiate the LocalFileStore with the temporary directory as the root path
        store = LocalFileStore(temp_dir)
        yield store


def test_create_lc_store(file_store: LocalFileStore) -> None:
    """Test that a docstore is created from a base store."""
    docstore = create_lc_store(file_store)
    docstore.mset([("key1", Document(page_content="hello", metadata={"key": "value"}))])
    fetched_doc = cast("Document", docstore.mget(["key1"])[0])
    assert fetched_doc.page_content == "hello"
    assert fetched_doc.metadata == {"key": "value"}


def test_create_kv_store(file_store: LocalFileStore) -> None:
    """Test that a docstore is created from a base store."""
    docstore = create_kv_docstore(file_store)
    docstore.mset([("key1", Document(page_content="hello", metadata={"key": "value"}))])
    fetched_doc = docstore.mget(["key1"])[0]
    assert isinstance(fetched_doc, Document)
    assert fetched_doc.page_content == "hello"
    assert fetched_doc.metadata == {"key": "value"}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
