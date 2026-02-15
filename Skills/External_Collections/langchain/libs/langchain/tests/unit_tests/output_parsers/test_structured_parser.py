# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Any

from langchain_core.exceptions import OutputParserException

from langchain_classic.output_parsers import ResponseSchema, StructuredOutputParser


def test_parse() -> None:
    """Test parsing structured output."""
    response_schemas = [
        ResponseSchema(name="name", description="desc"),
        ResponseSchema(name="age", description="desc"),
    ]
    parser = StructuredOutputParser.from_response_schemas(response_schemas)

    # Test valid JSON input
    text = '```json\n{"name": "John", "age": 30}\n```'
    expected_result = {"name": "John", "age": 30}
    result = parser.parse(text)
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

    # Test invalid JSON input
    text = '```json\n{"name": "John"}\n```'
    try:
        parser.parse(text)
    except OutputParserException:
        pass  # Test passes if OutputParserException is raised
    else:
        msg = f"Expected OutputParserException, but got {parser.parse(text)}"
        raise AssertionError(msg)


def test_output_type() -> None:
    """Test the output type of the structured output parser is Dict[str, Any]."""
    response_schemas = [
        ResponseSchema(name="name", description="desc"),
        ResponseSchema(name="age", description="desc"),
    ]
    parser = StructuredOutputParser.from_response_schemas(response_schemas)
    assert parser.OutputType == dict[str, Any]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
