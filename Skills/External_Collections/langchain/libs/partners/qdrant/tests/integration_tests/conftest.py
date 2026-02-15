# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os

from qdrant_client import QdrantClient

from tests.integration_tests.fixtures import qdrant_locations


def pytest_runtest_teardown() -> None:
    """Clean up all collections after the each test."""
    for location in qdrant_locations():
        client = QdrantClient(location=location, api_key=os.getenv("QDRANT_API_KEY"))
        collections = client.get_collections().collections

        for collection in collections:
            client.delete_collection(collection.name)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
