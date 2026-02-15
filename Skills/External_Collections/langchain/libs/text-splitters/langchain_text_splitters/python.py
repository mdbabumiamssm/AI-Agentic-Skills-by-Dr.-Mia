# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Python code text splitter."""

from __future__ import annotations

from typing import Any

from langchain_text_splitters.base import Language
from langchain_text_splitters.character import RecursiveCharacterTextSplitter


class PythonCodeTextSplitter(RecursiveCharacterTextSplitter):
    """Attempts to split the text along Python syntax."""

    def __init__(self, **kwargs: Any) -> None:
        """Initialize a PythonCodeTextSplitter."""
        separators = self.get_separators_for_language(Language.PYTHON)
        super().__init__(separators=separators, **kwargs)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
