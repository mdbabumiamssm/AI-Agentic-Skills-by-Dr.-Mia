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
import uuid

import pytest  # type: ignore[import-not-found]

from langchain_qdrant import Qdrant
from tests.integration_tests.common import ConsistentFakeEmbeddings


@pytest.mark.parametrize("vector_name", ["custom-vector"])
def test_qdrant_from_existing_collection_uses_same_collection(vector_name: str) -> None:
    """Test if the Qdrant.from_existing_collection reuses the same collection."""
    from qdrant_client import QdrantClient

    collection_name = uuid.uuid4().hex
    with tempfile.TemporaryDirectory() as tmpdir:
        docs = ["foo"]
        qdrant = Qdrant.from_texts(
            docs,
            embedding=ConsistentFakeEmbeddings(),
            path=str(tmpdir),
            collection_name=collection_name,
            vector_name=vector_name,
        )
        del qdrant

        qdrant = Qdrant.from_existing_collection(
            embedding=ConsistentFakeEmbeddings(),
            path=str(tmpdir),
            collection_name=collection_name,
            vector_name=vector_name,
        )
        qdrant.add_texts(["baz", "bar"])
        del qdrant

        client = QdrantClient(path=str(tmpdir))
        assert client.count(collection_name).count == 3

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
