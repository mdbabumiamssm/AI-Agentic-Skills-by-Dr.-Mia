# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from agbenchmark.utils.dependencies.graphs import is_circular


def test_is_circular():
    cyclic_graph = {
        "nodes": [
            {"id": "A", "data": {"category": []}},
            {"id": "B", "data": {"category": []}},
            {"id": "C", "data": {"category": []}},
            {"id": "D", "data": {"category": []}},  # New node
        ],
        "edges": [
            {"from": "A", "to": "B"},
            {"from": "B", "to": "C"},
            {"from": "C", "to": "D"},
            {"from": "D", "to": "A"},  # This edge creates a cycle
        ],
    }

    result = is_circular(cyclic_graph)
    assert result is not None, "Expected a cycle, but none was detected"
    assert all(
        (
            (result[i], result[i + 1])
            in [(x["from"], x["to"]) for x in cyclic_graph["edges"]]
        )
        for i in range(len(result) - 1)
    ), "The detected cycle path is not part of the graph's edges"


def test_is_not_circular():
    acyclic_graph = {
        "nodes": [
            {"id": "A", "data": {"category": []}},
            {"id": "B", "data": {"category": []}},
            {"id": "C", "data": {"category": []}},
            {"id": "D", "data": {"category": []}},  # New node
        ],
        "edges": [
            {"from": "A", "to": "B"},
            {"from": "B", "to": "C"},
            {"from": "C", "to": "D"},
            # No back edge from D to any node, so it remains acyclic
        ],
    }

    assert is_circular(acyclic_graph) is None, "Detected a cycle in an acyclic graph"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
