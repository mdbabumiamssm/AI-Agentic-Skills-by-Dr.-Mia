# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Integration tests for Exa find similar tool."""

from langchain_exa import (
    ExaFindSimilarResults,  # type: ignore[import-not-found, import-not-found]
)


def test_similarity_tool() -> None:
    """Test that the Exa find similar tool works."""
    tool = ExaFindSimilarResults()
    res = tool.invoke(
        {
            "url": "https://boutiquejapan.com/when-is-the-best-time-of-year-to-visit-japan/",
            "num_results": 5,
        }
    )
    print(res)  # noqa: T201
    assert not isinstance(res, str)  # str means error for this tool

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
