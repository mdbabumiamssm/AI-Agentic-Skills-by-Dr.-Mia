# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test Nomic embeddings."""

from langchain_nomic.embeddings import NomicEmbeddings


def test_langchain_nomic_embedding_documents() -> None:
    """Test nomic embeddings."""
    documents = ["foo bar"]
    embedding = NomicEmbeddings(model="nomic-embed-text-v1")
    output = embedding.embed_documents(documents)
    assert len(output) == 1
    assert len(output[0]) > 0


def test_langchain_nomic_embedding_query() -> None:
    """Test nomic embeddings."""
    document = "foo bar"
    embedding = NomicEmbeddings(model="nomic-embed-text-v1")
    output = embedding.embed_query(document)
    assert len(output) > 0


def test_langchain_nomic_embedding_dimensionality() -> None:
    """Test nomic embeddings."""
    documents = ["foo bar"]
    embedding = NomicEmbeddings(model="nomic-embed-text-v1.5", dimensionality=256)
    output = embedding.embed_documents(documents)
    assert len(output) == 1
    assert len(output[0]) == 256

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
