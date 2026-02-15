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
    from langchain_community.document_transformers import (
        EmbeddingsClusteringFilter,
        EmbeddingsRedundantFilter,
        get_stateful_documents,
    )
    from langchain_community.document_transformers.embeddings_redundant_filter import (
        _DocumentWithState,
        _filter_similar_embeddings,
        _get_embeddings_from_stateful_docs,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "EmbeddingsRedundantFilter": "langchain_community.document_transformers",
    "EmbeddingsClusteringFilter": "langchain_community.document_transformers",
    "_DocumentWithState": (
        "langchain_community.document_transformers.embeddings_redundant_filter"
    ),
    "get_stateful_documents": "langchain_community.document_transformers",
    "_get_embeddings_from_stateful_docs": (
        "langchain_community.document_transformers.embeddings_redundant_filter"
    ),
    "_filter_similar_embeddings": (
        "langchain_community.document_transformers.embeddings_redundant_filter"
    ),
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "EmbeddingsClusteringFilter",
    "EmbeddingsRedundantFilter",
    "_DocumentWithState",
    "_filter_similar_embeddings",
    "_get_embeddings_from_stateful_docs",
    "get_stateful_documents",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
