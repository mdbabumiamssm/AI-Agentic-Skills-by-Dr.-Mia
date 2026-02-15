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

import pytest
from langchain_core.exceptions import OutputParserException

from langchain_classic.output_parsers.enum import EnumOutputParser


class Colors(Enum):
    RED = "red"
    GREEN = "green"
    BLUE = "blue"


def test_enum_output_parser_parse() -> None:
    parser = EnumOutputParser(enum=Colors)

    # Test valid inputs
    result = parser.parse("red")
    assert result == Colors.RED

    result = parser.parse("green")
    assert result == Colors.GREEN

    result = parser.parse("blue")
    assert result == Colors.BLUE

    # Test invalid input
    with pytest.raises(OutputParserException):
        parser.parse("INVALID")


def test_enum_output_parser_output_type() -> None:
    """Test the output type of the enum output parser is the expected enum."""
    assert EnumOutputParser(enum=Colors).OutputType is Colors

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
