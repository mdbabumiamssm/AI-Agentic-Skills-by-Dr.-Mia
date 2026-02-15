# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import re

import pytest

from langchain_classic.output_parsers.boolean import BooleanOutputParser


def test_boolean_output_parser_parse() -> None:
    parser = BooleanOutputParser()

    # Test valid input
    result = parser.parse("YES")
    assert result is True

    # Test valid input
    result = parser.parse("NO")
    assert result is False

    # Test valid input
    result = parser.parse("yes")
    assert result is True

    # Test valid input
    result = parser.parse("no")
    assert result is False

    # Test valid input
    result = parser.parse("Not relevant (NO)")
    assert result is False

    # Test valid input
    result = parser.parse("NOW this is relevant (YES)")
    assert result is True

    # Test ambiguous input
    with pytest.raises(
        ValueError,
        match=re.escape("Ambiguous response. Both YES and NO in received: YES NO."),
    ):
        parser.parse("YES NO")

    with pytest.raises(
        ValueError,
        match=re.escape("Ambiguous response. Both YES and NO in received: NO YES."),
    ):
        parser.parse("NO YES")
    # Bad input
    with pytest.raises(
        ValueError,
        match=re.escape(
            "BooleanOutputParser expected output value to include either YES or NO. "
            "Received BOOM."
        ),
    ):
        parser.parse("BOOM")


def test_boolean_output_parser_output_type() -> None:
    """Test the output type of the boolean output parser is a boolean."""
    assert BooleanOutputParser().OutputType is bool

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
