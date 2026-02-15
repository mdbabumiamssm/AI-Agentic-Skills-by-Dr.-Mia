# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""**Graphs** provide a natural language interface to graph databases."""

from typing import TYPE_CHECKING, Any

from langchain_classic._api import create_importer

if TYPE_CHECKING:
    from langchain_community.graphs import (
        ArangoGraph,
        FalkorDBGraph,
        HugeGraph,
        KuzuGraph,
        MemgraphGraph,
        NebulaGraph,
        Neo4jGraph,
        NeptuneGraph,
        NetworkxEntityGraph,
        RdfGraph,
    )


# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "MemgraphGraph": "langchain_community.graphs",
    "NetworkxEntityGraph": "langchain_community.graphs",
    "Neo4jGraph": "langchain_community.graphs",
    "NebulaGraph": "langchain_community.graphs",
    "NeptuneGraph": "langchain_community.graphs",
    "KuzuGraph": "langchain_community.graphs",
    "HugeGraph": "langchain_community.graphs",
    "RdfGraph": "langchain_community.graphs",
    "ArangoGraph": "langchain_community.graphs",
    "FalkorDBGraph": "langchain_community.graphs",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "ArangoGraph",
    "FalkorDBGraph",
    "HugeGraph",
    "KuzuGraph",
    "MemgraphGraph",
    "NebulaGraph",
    "Neo4jGraph",
    "NeptuneGraph",
    "NetworkxEntityGraph",
    "RdfGraph",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
