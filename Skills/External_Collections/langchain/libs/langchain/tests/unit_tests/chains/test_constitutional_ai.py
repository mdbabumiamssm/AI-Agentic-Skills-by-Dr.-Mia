# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Unit tests for the Constitutional AI chain."""

from langchain_classic.chains.constitutional_ai.base import ConstitutionalChain

TEXT_ONE = """ This text is bad.

Revision request: Make it better.

Revision:"""

TEXT_TWO = """ This text is bad.\n\n"""

TEXT_THREE = """ This text is bad.

Revision request: Make it better.

Revision: Better text"""


def test_critique_parsing() -> None:
    """Test parsing of critique text."""
    for text in [TEXT_ONE, TEXT_TWO, TEXT_THREE]:
        critique = ConstitutionalChain._parse_critique(text)

        assert critique.strip() == "This text is bad.", (
            f"Failed on {text} with {critique}"
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
