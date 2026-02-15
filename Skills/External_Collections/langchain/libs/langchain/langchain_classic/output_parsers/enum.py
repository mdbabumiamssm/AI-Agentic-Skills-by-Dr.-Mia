# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from enum import Enum

from langchain_core.exceptions import OutputParserException
from langchain_core.output_parsers import BaseOutputParser
from langchain_core.utils import pre_init
from typing_extensions import override


class EnumOutputParser(BaseOutputParser[Enum]):
    """Parse an output that is one of a set of values."""

    enum: type[Enum]
    """The enum to parse. Its values must be strings."""

    @pre_init
    def _raise_deprecation(cls, values: dict) -> dict:
        enum = values["enum"]
        if not all(isinstance(e.value, str) for e in enum):
            msg = "Enum values must be strings"
            raise ValueError(msg)
        return values

    @property
    def _valid_values(self) -> list[str]:
        return [e.value for e in self.enum]

    @override
    def parse(self, response: str) -> Enum:
        try:
            return self.enum(response.strip())
        except ValueError as e:
            msg = (
                f"Response '{response}' is not one of the "
                f"expected values: {self._valid_values}"
            )
            raise OutputParserException(msg) from e

    @override
    def get_format_instructions(self) -> str:
        return f"Select one of the following options: {', '.join(self._valid_values)}"

    @property
    @override
    def OutputType(self) -> type[Enum]:
        return self.enum

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
