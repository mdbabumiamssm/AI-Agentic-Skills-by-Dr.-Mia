# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic._api.module_import import create_importer


def test_import_from_non_deprecated_path() -> None:
    """Test importing all modules in langchain."""
    module_lookup = {
        "Document": "langchain_core.documents",
    }
    lookup = create_importer(__package__, module_lookup=module_lookup)
    imported_doc = lookup("Document")
    from langchain_core.documents import Document

    assert imported_doc is Document


def test_import_from_deprecated_path() -> None:
    """Test importing all modules in langchain."""
    module_lookup = {
        "Document": "langchain_core.documents",
    }
    lookup = create_importer(__package__, deprecated_lookups=module_lookup)
    imported_doc = lookup("Document")

    from langchain_core.documents import Document

    assert imported_doc is Document


def test_import_using_fallback_module() -> None:
    """Test import using fallback module."""
    lookup = create_importer(__package__, fallback_module="langchain_core.documents")
    imported_doc = lookup("Document")
    from langchain_core.documents import Document

    assert imported_doc is Document

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
