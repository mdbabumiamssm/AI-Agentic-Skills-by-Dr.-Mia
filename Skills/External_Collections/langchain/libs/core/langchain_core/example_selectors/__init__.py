# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Example selectors.

**Example selector** implements logic for selecting examples to include them in prompts.
This allows us to select examples that are most relevant to the input.
"""

from typing import TYPE_CHECKING

from langchain_core._import_utils import import_attr

if TYPE_CHECKING:
    from langchain_core.example_selectors.base import BaseExampleSelector
    from langchain_core.example_selectors.length_based import (
        LengthBasedExampleSelector,
    )
    from langchain_core.example_selectors.semantic_similarity import (
        MaxMarginalRelevanceExampleSelector,
        SemanticSimilarityExampleSelector,
        sorted_values,
    )

__all__ = (
    "BaseExampleSelector",
    "LengthBasedExampleSelector",
    "MaxMarginalRelevanceExampleSelector",
    "SemanticSimilarityExampleSelector",
    "sorted_values",
)

_dynamic_imports = {
    "BaseExampleSelector": "base",
    "LengthBasedExampleSelector": "length_based",
    "MaxMarginalRelevanceExampleSelector": "semantic_similarity",
    "SemanticSimilarityExampleSelector": "semantic_similarity",
    "sorted_values": "semantic_similarity",
}


def __getattr__(attr_name: str) -> object:
    module_name = _dynamic_imports.get(attr_name)
    result = import_attr(attr_name, module_name, __spec__.parent)
    globals()[attr_name] = result
    return result


def __dir__() -> list[str]:
    return list(__all__)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
