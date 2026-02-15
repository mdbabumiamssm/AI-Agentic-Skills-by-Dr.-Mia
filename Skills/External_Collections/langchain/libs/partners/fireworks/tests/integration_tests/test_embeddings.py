# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test Fireworks embeddings."""

from langchain_fireworks.embeddings import FireworksEmbeddings


def test_langchain_fireworks_embedding_documents() -> None:
    """Test Fireworks hosted embeddings."""
    documents = ["foo bar"]
    embedding = FireworksEmbeddings(model="nomic-ai/nomic-embed-text-v1.5")
    output = embedding.embed_documents(documents)
    assert len(output) == 1
    assert len(output[0]) > 0


def test_langchain_fireworks_embedding_query() -> None:
    """Test Fireworks hosted embeddings."""
    document = "foo bar"
    embedding = FireworksEmbeddings(model="nomic-ai/nomic-embed-text-v1.5")
    output = embedding.embed_query(document)
    assert len(output) > 0

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
