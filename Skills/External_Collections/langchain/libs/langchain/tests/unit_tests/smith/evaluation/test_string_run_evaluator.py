# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Tests for the string run evaluator."""

from unittest.mock import MagicMock

from langchain_classic.evaluation import criteria
from langchain_classic.smith.evaluation.string_run_evaluator import (
    ChainStringRunMapper,
    StringRunEvaluatorChain,
)
from tests.unit_tests.llms import fake_llm


def test_evaluate_run() -> None:
    run_mapper = ChainStringRunMapper()
    string_evaluator = criteria.CriteriaEvalChain.from_llm(fake_llm.FakeLLM())
    evaluator = StringRunEvaluatorChain(
        run_mapper=run_mapper,
        example_mapper=None,
        name="test_evaluator",
        string_evaluator=string_evaluator,
    )
    run = MagicMock()
    example = MagicMock()
    res = evaluator.evaluate_run(run, example)
    assert str(res.comment).startswith("Error evaluating run ")
    assert res.key == string_evaluator.evaluation_name

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
