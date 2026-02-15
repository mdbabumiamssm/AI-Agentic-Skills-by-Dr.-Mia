# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.indexing import __all__


def test_all() -> None:
    """Use to catch obvious breaking changes."""
    assert list(__all__) == sorted(__all__, key=str)
    assert set(__all__) == {
        "aindex",
        "DeleteResponse",
        "DocumentIndex",
        "index",
        "IndexingResult",
        "InMemoryRecordManager",
        "RecordManager",
        "UpsertResponse",
    }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
