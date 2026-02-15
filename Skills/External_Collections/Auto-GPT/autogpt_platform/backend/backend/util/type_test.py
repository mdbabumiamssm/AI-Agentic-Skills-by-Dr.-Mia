# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import List, Optional

from backend.util.type import convert


def test_type_conversion():
    assert convert(5.5, int) == 5
    assert convert("5.5", int) == 5
    assert convert([1, 2, 3], int) == 3
    assert convert("7", Optional[int]) == 7
    assert convert("7", int | None) == 7

    assert convert("5.5", float) == 5.5
    assert convert(5, float) == 5.0

    assert convert("True", bool) is True
    assert convert("False", bool) is False

    assert convert(5, str) == "5"
    assert convert({"a": 1, "b": 2}, str) == '{"a": 1, "b": 2}'
    assert convert([1, 2, 3], str) == "[1, 2, 3]"

    assert convert("5", list) == ["5"]
    assert convert((1, 2, 3), list) == [1, 2, 3]
    assert convert({1, 2, 3}, list) == [1, 2, 3]

    assert convert("5", dict) == {"value": 5}
    assert convert('{"a": 1, "b": 2}', dict) == {"a": 1, "b": 2}
    assert convert([1, 2, 3], dict) == {0: 1, 1: 2, 2: 3}
    assert convert((1, 2, 3), dict) == {0: 1, 1: 2, 2: 3}

    assert convert("5", List[int]) == [5]
    assert convert("[5,4,2]", List[int]) == [5, 4, 2]
    assert convert([5, 4, 2], List[str]) == ["5", "4", "2"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
