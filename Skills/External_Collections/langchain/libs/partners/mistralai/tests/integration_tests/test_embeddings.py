# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test MistralAI Embedding."""

from unittest.mock import patch

import httpx
import pytest
import tenacity

from langchain_mistralai import MistralAIEmbeddings


def test_mistralai_embedding_documents() -> None:
    """Test MistralAI embeddings for documents."""
    documents = ["foo bar", "test document"]
    embedding = MistralAIEmbeddings()
    output = embedding.embed_documents(documents)
    assert len(output) == 2
    assert len(output[0]) == 1024


def test_mistralai_embedding_query() -> None:
    """Test MistralAI embeddings for query."""
    document = "foo bar"
    embedding = MistralAIEmbeddings()
    output = embedding.embed_query(document)
    assert len(output) == 1024


async def test_mistralai_embedding_documents_async() -> None:
    """Test MistralAI embeddings for documents."""
    documents = ["foo bar", "test document"]
    embedding = MistralAIEmbeddings()
    output = await embedding.aembed_documents(documents)
    assert len(output) == 2
    assert len(output[0]) == 1024


async def test_mistralai_embedding_documents_tenacity_error_async() -> None:
    """Test MistralAI embeddings for documents."""
    documents = ["foo bar", "test document"]
    embedding = MistralAIEmbeddings(max_retries=0)
    mock_response = httpx.Response(
        status_code=400,
        request=httpx.Request("POST", url=embedding.async_client.base_url),
    )
    with (
        patch.object(embedding.async_client, "post", return_value=mock_response),
        pytest.raises(tenacity.RetryError),
    ):
        await embedding.aembed_documents(documents)


async def test_mistralai_embedding_documents_http_error_async() -> None:
    """Test MistralAI embeddings for documents."""
    documents = ["foo bar", "test document"]
    embedding = MistralAIEmbeddings(max_retries=None)
    mock_response = httpx.Response(
        status_code=400,
        request=httpx.Request("POST", url=embedding.async_client.base_url),
    )
    with (
        patch.object(embedding.async_client, "post", return_value=mock_response),
        pytest.raises(httpx.HTTPStatusError),
    ):
        await embedding.aembed_documents(documents)


async def test_mistralai_embedding_query_async() -> None:
    """Test MistralAI embeddings for query."""
    document = "foo bar"
    embedding = MistralAIEmbeddings()
    output = await embedding.aembed_query(document)
    assert len(output) == 1024


def test_mistralai_embedding_documents_long() -> None:
    """Test MistralAI embeddings for documents."""
    documents = ["foo bar " * 1000, "test document " * 1000] * 5
    embedding = MistralAIEmbeddings()
    output = embedding.embed_documents(documents)
    assert len(output) == 10
    assert len(output[0]) == 1024


def test_mistralai_embed_query_character() -> None:
    """Test MistralAI embeddings for query."""
    document = "😳"
    embedding = MistralAIEmbeddings()
    output = embedding.embed_query(document)
    assert len(output) == 1024

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
