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
    from langchain_community.document_loaders import (
        UnstructuredAPIFileIOLoader,
        UnstructuredAPIFileLoader,
        UnstructuredFileIOLoader,
        UnstructuredFileLoader,
    )
    from langchain_community.document_loaders.unstructured import (
        UnstructuredBaseLoader,
        get_elements_from_api,
        satisfies_min_unstructured_version,
        validate_unstructured_version,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "satisfies_min_unstructured_version": (
        "langchain_community.document_loaders.unstructured"
    ),
    "validate_unstructured_version": (
        "langchain_community.document_loaders.unstructured"
    ),
    "UnstructuredBaseLoader": "langchain_community.document_loaders.unstructured",
    "UnstructuredFileLoader": "langchain_community.document_loaders",
    "get_elements_from_api": "langchain_community.document_loaders.unstructured",
    "UnstructuredAPIFileLoader": "langchain_community.document_loaders",
    "UnstructuredFileIOLoader": "langchain_community.document_loaders",
    "UnstructuredAPIFileIOLoader": "langchain_community.document_loaders",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "UnstructuredAPIFileIOLoader",
    "UnstructuredAPIFileLoader",
    "UnstructuredBaseLoader",
    "UnstructuredFileIOLoader",
    "UnstructuredFileLoader",
    "get_elements_from_api",
    "satisfies_min_unstructured_version",
    "validate_unstructured_version",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
