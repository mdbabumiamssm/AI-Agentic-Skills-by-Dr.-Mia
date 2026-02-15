# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

r"""Comparison evaluators.

This module contains evaluators for comparing the output of two models,
be they LLMs, Chains, or otherwise. This can be used for scoring
preferences, measuring similarity / semantic equivalence between outputs,
or any other comparison task.

Example:
    >>> from langchain_openai import ChatOpenAI
    >>> from langchain_classic.evaluation.comparison import PairwiseStringEvalChain
    >>> llm = ChatOpenAI(temperature=0)
    >>> chain = PairwiseStringEvalChain.from_llm(llm=llm)
    >>> result = chain.evaluate_string_pairs(
    ...     input = "What is the chemical formula for water?",
    ...     prediction = "H2O",
    ...     prediction_b = (
    ...        "The chemical formula for water is H2O, which means"
    ...        " there are two hydrogen atoms and one oxygen atom."
    ...     reference = "The chemical formula for water is H2O.",
    ... )
    >>> print(result)
    # {
    #    "value": "B",
    #    "comment": "Both responses accurately state"
    #       " that the chemical formula for water is H2O."
    #       " However, Response B provides additional information"
    # .     " by explaining what the formula means.\n[[B]]"
    # }
"""

from langchain_classic.evaluation.comparison.eval_chain import (
    LabeledPairwiseStringEvalChain,
    PairwiseStringEvalChain,
)

__all__ = ["LabeledPairwiseStringEvalChain", "PairwiseStringEvalChain"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
