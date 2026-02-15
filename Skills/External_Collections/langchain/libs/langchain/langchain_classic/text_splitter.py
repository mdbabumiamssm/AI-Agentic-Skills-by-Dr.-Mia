# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Kept for backwards compatibility."""

from langchain_text_splitters import (
    Language,
    RecursiveCharacterTextSplitter,
    TextSplitter,
    Tokenizer,
    TokenTextSplitter,
)
from langchain_text_splitters.base import split_text_on_tokens
from langchain_text_splitters.character import CharacterTextSplitter
from langchain_text_splitters.html import ElementType, HTMLHeaderTextSplitter
from langchain_text_splitters.json import RecursiveJsonSplitter
from langchain_text_splitters.konlpy import KonlpyTextSplitter
from langchain_text_splitters.latex import LatexTextSplitter
from langchain_text_splitters.markdown import (
    HeaderType,
    LineType,
    MarkdownHeaderTextSplitter,
    MarkdownTextSplitter,
)
from langchain_text_splitters.nltk import NLTKTextSplitter
from langchain_text_splitters.python import PythonCodeTextSplitter
from langchain_text_splitters.sentence_transformers import (
    SentenceTransformersTokenTextSplitter,
)
from langchain_text_splitters.spacy import SpacyTextSplitter

__all__ = [
    "CharacterTextSplitter",
    "ElementType",
    "HTMLHeaderTextSplitter",
    "HeaderType",
    "KonlpyTextSplitter",
    "Language",
    "LatexTextSplitter",
    "LineType",
    "MarkdownHeaderTextSplitter",
    "MarkdownTextSplitter",
    "NLTKTextSplitter",
    "PythonCodeTextSplitter",
    "RecursiveCharacterTextSplitter",
    "RecursiveJsonSplitter",
    "SentenceTransformersTokenTextSplitter",
    "SpacyTextSplitter",
    "TextSplitter",
    "TokenTextSplitter",
    "Tokenizer",
    "split_text_on_tokens",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
