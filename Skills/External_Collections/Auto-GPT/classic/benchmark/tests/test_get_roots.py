# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from agbenchmark.utils.dependencies.graphs import get_roots


def test_get_roots():
    graph = {
        "nodes": [
            {"id": "A", "data": {"category": []}},
            {"id": "B", "data": {"category": []}},
            {"id": "C", "data": {"category": []}},
            {"id": "D", "data": {"category": []}},
        ],
        "edges": [
            {"from": "A", "to": "B"},
            {"from": "B", "to": "C"},
        ],
    }

    result = get_roots(graph)
    assert set(result) == {
        "A",
        "D",
    }, f"Expected roots to be 'A' and 'D', but got {result}"


def test_no_roots():
    fully_connected_graph = {
        "nodes": [
            {"id": "A", "data": {"category": []}},
            {"id": "B", "data": {"category": []}},
            {"id": "C", "data": {"category": []}},
        ],
        "edges": [
            {"from": "A", "to": "B"},
            {"from": "B", "to": "C"},
            {"from": "C", "to": "A"},
        ],
    }

    result = get_roots(fully_connected_graph)
    assert not result, "Expected no roots, but found some"


# def test_no_rcoots():
#     fully_connected_graph = {
#         "nodes": [
#             {"id": "A", "data": {"category": []}},
#             {"id": "B", "data": {"category": []}},
#             {"id": "C", "data": {"category": []}},
#         ],
#         "edges": [
#             {"from": "A", "to": "B"},
#             {"from": "D", "to": "C"},
#         ],
#     }
#
#     result = get_roots(fully_connected_graph)
#     assert set(result) == {"A"}, f"Expected roots to be 'A', but got {result}"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
