# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic.document_loaders.parsers import __all__


def test_parsers_public_api_correct() -> None:
    """Test public API of parsers for breaking changes."""
    assert set(__all__) == {
        "BS4HTMLParser",
        "DocAIParser",
        "GrobidParser",
        "LanguageParser",
        "OpenAIWhisperParser",
        "PyPDFParser",
        "PDFMinerParser",
        "PyMuPDFParser",
        "PyPDFium2Parser",
        "PDFPlumberParser",
    }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
