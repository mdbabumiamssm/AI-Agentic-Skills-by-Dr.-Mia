# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest

from langchain_classic.evaluation import ExactMatchStringEvaluator


@pytest.fixture
def exact_match_string_evaluator() -> ExactMatchStringEvaluator:
    """Create an ExactMatchStringEvaluator with default configuration."""
    return ExactMatchStringEvaluator()


@pytest.fixture
def exact_match_string_evaluator_ignore_case() -> ExactMatchStringEvaluator:
    """Create an ExactMatchStringEvaluator with ignore_case set to True."""
    return ExactMatchStringEvaluator(ignore_case=True)


def test_default_exact_matching(
    exact_match_string_evaluator: ExactMatchStringEvaluator,
) -> None:
    prediction = "Mindy is the CTO"
    reference = "Mindy is the CTO"
    result = exact_match_string_evaluator.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 1.0

    reference = "Mindy is the CEO"
    result = exact_match_string_evaluator.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 0.0


def test_exact_matching_with_ignore_case(
    exact_match_string_evaluator_ignore_case: ExactMatchStringEvaluator,
) -> None:
    prediction = "Mindy is the CTO"
    reference = "mindy is the cto"
    result = exact_match_string_evaluator_ignore_case.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 1.0

    reference = "mindy is the CEO"
    result = exact_match_string_evaluator_ignore_case.evaluate_strings(
        prediction=prediction,
        reference=reference,
    )
    assert result["score"] == 0.0

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
