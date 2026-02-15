# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Unit tests for imports in `langchain_exa`."""

from langchain_exa import __all__  # type: ignore[import-not-found, import-not-found]

EXPECTED_ALL = [
    "ExaSearchResults",
    "ExaSearchRetriever",
    "HighlightsContentsOptions",
    "TextContentsOptions",
    "ExaFindSimilarResults",
]


def test_all_imports() -> None:
    """Test that all expected imports are in `__all__`."""
    assert sorted(EXPECTED_ALL) == sorted(__all__)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
