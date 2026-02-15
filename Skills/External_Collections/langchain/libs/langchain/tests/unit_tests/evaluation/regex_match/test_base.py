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

from langchain_classic.evaluation import RegexMatchStringEvaluator


@pytest.fixture
def regex_match_string_evaluator() -> RegexMatchStringEvaluator:
    """Create a RegexMatchStringEvaluator with default configuration."""
    return RegexMatchStringEvaluator()


@pytest.fixture
def regex_match_string_evaluator_ignore_case() -> RegexMatchStringEvaluator:
    """Create a RegexMatchStringEvaluator with IGNORECASE flag."""
    return RegexMatchStringEvaluator(flags=re.IGNORECASE)


def test_default_regex_matching(
    regex_match_string_evaluator: RegexMatchStringEvaluator,
) -> None:
    prediction = "Mindy is the CTO"
    reference = "^Mindy.*CTO$"
    result = regex_match_string_evaluator.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 1.0

    reference = "^Mike.*CEO$"
    result = regex_match_string_evaluator.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 0.0


def test_regex_matching_with_ignore_case(
    regex_match_string_evaluator_ignore_case: RegexMatchStringEvaluator,
) -> None:
    prediction = "Mindy is the CTO"
    reference = "^mindy.*cto$"
    result = regex_match_string_evaluator_ignore_case.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 1.0

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
