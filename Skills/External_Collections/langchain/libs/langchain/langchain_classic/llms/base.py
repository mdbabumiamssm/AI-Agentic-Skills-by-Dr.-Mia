# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""This module provides backward-compatible exports of core language model classes.

These classes are re-exported for compatibility with older versions of LangChain
and allow users to import language model interfaces from a stable path.

Exports:
    - LLM: Abstract base class for all LLMs
    - BaseLLM: Deprecated or foundational class for legacy LLMs
    - BaseLanguageModel: Base class for core language model implementations
"""

from langchain_core.language_models import BaseLanguageModel
from langchain_core.language_models.llms import LLM, BaseLLM

__all__ = [
    "LLM",
    "BaseLLM",
    "BaseLanguageModel",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
