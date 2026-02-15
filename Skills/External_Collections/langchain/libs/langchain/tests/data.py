# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Module defines common test data."""

from pathlib import Path

_THIS_DIR = Path(__file__).parent

_EXAMPLES_DIR = _THIS_DIR / "integration_tests" / "examples"

# Paths to test PDF files
HELLO_PDF = _EXAMPLES_DIR / "hello.pdf"
LAYOUT_PARSER_PAPER_PDF = _EXAMPLES_DIR / "layout-parser-paper.pdf"
DUPLICATE_CHARS = _EXAMPLES_DIR / "duplicate-chars.pdf"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
