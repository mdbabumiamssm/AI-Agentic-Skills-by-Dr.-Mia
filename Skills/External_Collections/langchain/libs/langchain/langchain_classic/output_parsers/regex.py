# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

import re

from langchain_core.output_parsers import BaseOutputParser
from typing_extensions import override


class RegexParser(BaseOutputParser[dict[str, str]]):
    """Parse the output of an LLM call using a regex."""

    @classmethod
    @override
    def is_lc_serializable(cls) -> bool:
        return True

    regex: str
    """The regex to use to parse the output."""
    output_keys: list[str]
    """The keys to use for the output."""
    default_output_key: str | None = None
    """The default key to use for the output."""

    @property
    def _type(self) -> str:
        """Return the type key."""
        return "regex_parser"

    def parse(self, text: str) -> dict[str, str]:
        """Parse the output of an LLM call."""
        match = re.search(self.regex, text)
        if match:
            return {key: match.group(i + 1) for i, key in enumerate(self.output_keys)}
        if self.default_output_key is None:
            msg = f"Could not parse output: {text}"
            raise ValueError(msg)
        return {
            key: text if key == self.default_output_key else ""
            for key in self.output_keys
        }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
