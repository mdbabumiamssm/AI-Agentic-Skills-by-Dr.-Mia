# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Scoring evaluators.

This module contains evaluators for scoring on a 1-10 the output of models,
be they LLMs, Chains, or otherwise. This can be based on a variety of
criteria and or a reference answer.

Example:
    >>> from langchain_openai import ChatOpenAI
    >>> from langchain_classic.evaluation.scoring import ScoreStringEvalChain
    >>> model = ChatOpenAI(temperature=0, model_name="gpt-4")
    >>> chain = ScoreStringEvalChain.from_llm(llm=model)
    >>> result = chain.evaluate_strings(
    ...     input="What is the chemical formula for water?",
    ...     prediction="H2O",
    ...     reference="The chemical formula for water is H2O.",
    ... )
    >>> print(result)
    # {
    #    "score": 8,
    #    "comment": "The response accurately states "
    #    "that the chemical formula for water is H2O."
    #    "However, it does not provide an explanation of what the formula means."
    # }
"""

from langchain_classic.evaluation.scoring.eval_chain import (
    LabeledScoreStringEvalChain,
    ScoreStringEvalChain,
)

__all__ = ["LabeledScoreStringEvalChain", "ScoreStringEvalChain"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
