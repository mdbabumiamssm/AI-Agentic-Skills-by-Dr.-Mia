# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""LangChain integration for Exa."""

from exa_py.api import (
    HighlightsContentsOptions,
    TextContentsOptions,
)

from langchain_exa.retrievers import ExaSearchRetriever
from langchain_exa.tools import ExaFindSimilarResults, ExaSearchResults

__all__ = [
    "ExaFindSimilarResults",
    "ExaSearchResults",
    "ExaSearchRetriever",
    "HighlightsContentsOptions",
    "TextContentsOptions",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
