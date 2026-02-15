# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test map_rerank parser."""

import pytest

from langchain_classic.chains.question_answering.map_rerank_prompt import output_parser

GOOD_SCORE = "foo bar answer.\nScore: 80"
SCORE_WITH_EXPLANATION = (
    "foo bar answer.\n"
    "Score: 80 (fully answers the question, "
    "but could provide more detail on the specific error message)"
)


@pytest.mark.parametrize("answer", [GOOD_SCORE, SCORE_WITH_EXPLANATION])
def test_parse_scores(answer: str) -> None:
    result = output_parser.parse(answer)

    assert result["answer"] == "foo bar answer."

    score = int(result["score"])
    assert score == 80

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
