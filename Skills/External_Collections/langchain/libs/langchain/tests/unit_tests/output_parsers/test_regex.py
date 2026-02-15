# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic.output_parsers.regex import RegexParser

# NOTE: The almost same constant variables in ./test_combining_parser.py
DEF_EXPECTED_RESULT = {
    "confidence": "A",
    "explanation": "Paris is the capital of France according to Wikipedia.",
}

DEF_README = """```json
{
    "answer": "Paris",
    "source": "https://en.wikipedia.org/wiki/France"
}
```

//Confidence: A, Explanation: Paris is the capital of France according to Wikipedia."""


def test_regex_parser_parse() -> None:
    """Test regex parser parse."""
    parser = RegexParser(
        regex=r"Confidence: (A|B|C), Explanation: (.*)",
        output_keys=["confidence", "explanation"],
        default_output_key="noConfidence",
    )
    assert parser.parse(DEF_README) == DEF_EXPECTED_RESULT


def test_regex_parser_output_type() -> None:
    """Test regex parser output type is Dict[str, str]."""
    parser = RegexParser(
        regex=r"Confidence: (A|B|C), Explanation: (.*)",
        output_keys=["confidence", "explanation"],
        default_output_key="noConfidence",
    )
    assert parser.OutputType == dict[str, str]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
