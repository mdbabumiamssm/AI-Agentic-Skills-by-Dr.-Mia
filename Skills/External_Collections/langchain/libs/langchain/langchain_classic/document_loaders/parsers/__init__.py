# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import TYPE_CHECKING, Any

from langchain_classic._api import create_importer

if TYPE_CHECKING:
    from langchain_community.document_loaders.parsers.audio import OpenAIWhisperParser
    from langchain_community.document_loaders.parsers.docai import DocAIParser
    from langchain_community.document_loaders.parsers.grobid import GrobidParser
    from langchain_community.document_loaders.parsers.html.bs4 import BS4HTMLParser
    from langchain_community.document_loaders.parsers.language.language_parser import (
        LanguageParser,
    )
    from langchain_community.document_loaders.parsers.pdf import (
        PDFMinerParser,
        PDFPlumberParser,
        PyMuPDFParser,
        PyPDFium2Parser,
        PyPDFParser,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "BS4HTMLParser": "langchain_community.document_loaders.parsers.html.bs4",
    "DocAIParser": "langchain_community.document_loaders.parsers.docai",
    "GrobidParser": "langchain_community.document_loaders.parsers.grobid",
    "LanguageParser": (
        "langchain_community.document_loaders.parsers.language.language_parser"
    ),
    "OpenAIWhisperParser": "langchain_community.document_loaders.parsers.audio",
    "PDFMinerParser": "langchain_community.document_loaders.parsers.pdf",
    "PDFPlumberParser": "langchain_community.document_loaders.parsers.pdf",
    "PyMuPDFParser": "langchain_community.document_loaders.parsers.pdf",
    "PyPDFium2Parser": "langchain_community.document_loaders.parsers.pdf",
    "PyPDFParser": "langchain_community.document_loaders.parsers.pdf",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "BS4HTMLParser",
    "DocAIParser",
    "GrobidParser",
    "LanguageParser",
    "OpenAIWhisperParser",
    "PDFMinerParser",
    "PDFPlumberParser",
    "PyMuPDFParser",
    "PyPDFParser",
    "PyPDFium2Parser",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
