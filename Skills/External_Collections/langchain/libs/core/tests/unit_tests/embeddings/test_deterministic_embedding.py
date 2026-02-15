# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.embeddings import DeterministicFakeEmbedding


def test_deterministic_fake_embeddings() -> None:
    """Test that DeterministicFakeEmbedding is deterministic.

    Test that the deterministic fake embeddings return the same
    embedding vector for the same text.
    """
    fake = DeterministicFakeEmbedding(size=10)
    text = "Hello world!"
    assert fake.embed_query(text) == fake.embed_query(text)
    assert fake.embed_query(text) != fake.embed_query("Goodbye world!")
    assert fake.embed_documents([text, text]) == fake.embed_documents([text, text])
    assert fake.embed_documents([text, text]) != fake.embed_documents(
        [
            text,
            "Goodbye world!",
        ]
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
