# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""JavaScript framework text splitter."""

import re
from typing import Any

from langchain_text_splitters import RecursiveCharacterTextSplitter


class JSFrameworkTextSplitter(RecursiveCharacterTextSplitter):
    """Text splitter that handles React (JSX), Vue, and Svelte code.

    This splitter extends RecursiveCharacterTextSplitter to handle
    React (JSX), Vue, and Svelte code by:

    1. Detecting and extracting custom component tags from the text
    2. Using those tags as additional separators along with standard JS syntax

    The splitter combines:

    * Custom component tags as separators (e.g. <Component, <div)
    * JavaScript syntax elements (function, const, if, etc)
    * Standard text splitting on newlines

    This allows chunks to break at natural boundaries in
    React, Vue, and Svelte component code.
    """

    def __init__(
        self,
        separators: list[str] | None = None,
        chunk_size: int = 2000,
        chunk_overlap: int = 0,
        **kwargs: Any,
    ) -> None:
        """Initialize the JS Framework text splitter.

        Args:
            separators: Optional list of custom separator strings to use
            chunk_size: Maximum size of chunks to return
            chunk_overlap: Overlap in characters between chunks
            **kwargs: Additional arguments to pass to parent class
        """
        super().__init__(chunk_size=chunk_size, chunk_overlap=chunk_overlap, **kwargs)
        self._separators = separators or []

    def split_text(self, text: str) -> list[str]:
        """Split text into chunks.

        This method splits the text into chunks by:

        * Extracting unique opening component tags using regex
        * Creating separators list with extracted tags and JS separators
        * Splitting the text using the separators by calling the parent class method

        Args:
            text: String containing code to split

        Returns:
            List of text chunks split on component and JS boundaries
        """
        # Extract unique opening component tags using regex
        # Regex to match opening tags, excluding self-closing tags
        opening_tags = re.findall(r"<\s*([a-zA-Z0-9]+)[^>]*>", text)

        component_tags = []
        for tag in opening_tags:
            if tag not in component_tags:
                component_tags.append(tag)
        component_separators = [f"<{tag}" for tag in component_tags]

        js_separators = [
            "\nexport ",
            " export ",
            "\nfunction ",
            "\nasync function ",
            " async function ",
            "\nconst ",
            "\nlet ",
            "\nvar ",
            "\nclass ",
            " class ",
            "\nif ",
            " if ",
            "\nfor ",
            " for ",
            "\nwhile ",
            " while ",
            "\nswitch ",
            " switch ",
            "\ncase ",
            " case ",
            "\ndefault ",
            " default ",
        ]
        separators = (
            self._separators
            + js_separators
            + component_separators
            + ["<>", "\n\n", "&&\n", "||\n"]
        )
        self._separators = separators
        return super().split_text(text)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
